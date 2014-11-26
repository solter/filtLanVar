//============================================================================
// Routines for computing extreme or interior eigenvalues of a symmetric
// matrix by a filtered Lanczos procedure.
//
// Reference:
// "A Filtered Lanczos Procedure for Extreme and Interior Eigenvalue Problems"
// H.-r. Fang and Y. Saad, University of Minnesota Technical Report, 2011.
//============================================================================

#include <stdlib.h>  // for exit
#include <time.h>    // for time_t, time, clock_t, clock, CLOCKS_PER_SEC, etc.
#include <iostream>  // for cout, cerr, endl, etc. (under namespace std)

#include "matkit.h"

#include "laneig.h"
#include "polyfilt.h"
#include "filtlan.h"

using std::endl;


// this routine computes the eigenvalues in the requested window and the corresponding eigenvectors
// the window can contain the largest or smallest eigenvalues, in which case the requested are the extreme eigenvalues
// the window can also be well inside the range of the spectrum, in which case it is often referred to as an interior eigenvalue problem
// eigVal is a vector which contains the computed eigenvalues (if opts.wantEigVal==true)
// eigVec is a matrix whose columns are computed eigenvectors
//     more precisely, eigVec(:,j) is the eigenvector corresponding to the eigenvalue eigVal(j)
// S is the sparse matrix whose eigenvalues / eigenvectors are to be sought
// frame is a vector of 2 or 4 elements
//     if frame consists of 2 elements, then [frame(1),frame(2)] is the interval in which the eigenvalues are sought
//     if frame consists of 4 elements, then [frame(2),frame(3)] is the interval in which the eigenvalues are sought, and
//                                           [frame(1),frame(4)] is the interval which (tightly) contains all eigenvalues
//     in the former case, an interval which contains (tightly) all eigenvalues is not given, and will be determined by a standard Lanczos procedure
// polyDegree is the degree of s(z), with z*s(z) the polynomial filter
// baseDegree is the left-and-right degree of the base filter for each interval
// opts is a collection of options for the filtered Lanczos procedure
// the return FilteredLanczosInfo gives the information of the filtered Lanczos procedure
FilteredLanczosInfo FilteredLanczosEigenSolver(Vector &eigVal, Matrix &eigVec,
                                               const SparseMatrix &S, Vector &frame,
                                               mkIndex polyDegree, mkIndex baseDegree, FilteredLanczosOptions &opts) {
    // output information
    FilteredLanczosInfo outputInfo;

    // if the eigenvalue range is not given, determine it by the standard Lanczos process
    Real aa, bb;
    if (frame.Length() == 2) {
        // the window of desired eigenvalues is given as [frame(1), frame(2)]
        aa = frame(1);
        bb = frame(2);
        if (aa >= bb) {
            *opts.err << "FilteredLanczosEigenSolver(Vector&, Matrix&, const SparseMatrix&, Vector&, mkIndex, mkIndex, FilteredLanczosOptions&): ";
            *opts.err << "the requested interval of eigenvalues [" << aa << ", " << bb << "] should not have " << aa << ">=" << bb << "!" << endl;
            Basics::quit(1);
        }
        // an interval which (tightly) contains all eigenvalues is to be determined now
        clock_t rangeStart = clock();
        const char eigPart[] = "BE";
        LanczosOptions opts0;                       // default Lanczos options
        opts0.wantEigVec = true;                    // in case the default is false
        if (opts.numIterForEigenRange == 0u)
            opts.numIterForEigenRange = 20;         // set a default value
        opts0.minIter = opts.numIterForEigenRange;  // want exactly opts.numIterForEigenRange Lanczos iterations
        opts0.maxIter = opts.numIterForEigenRange;  // want exactly opts.numIterForEigenRange Lanczos iterations
        opts0.disp = 0;
        Vector eigVal0;
        Matrix eigVec0;
        outputInfo.forEigenRangeInfo = LanczosEigenSolver(eigVal0, eigVec0, S, 2, eigPart, opts0);
        Real lowerBound = eigVal0(1) - (S*eigVec0.column(1)-eigVal0(1)*eigVec0.column(1)).norm2();
        Real upperBound = eigVal0(2) + (S*eigVec0.column(2)-eigVal0(2)*eigVec0.column(2)).norm2();
        if (bb <= lowerBound || aa >= upperBound) {
            *opts.err << "FilteredLanczosEigenSolver(Vector&, Matrix&, const SparseMatrix&, Vector&, mkIndex, mkIndex, FilteredLanczosOptions&): ";
            *opts.err << "the requested interval of eigenvalues [" << aa << ", " << bb << "]" << " does not intersect the spectrum bounded in [" << lowerBound << ", " << upperBound << "]!" << endl;
            Basics::quit(1);
        }
        frame.resize(4);
        frame(1) = lowerBound;
        aa = Basics::max(lowerBound, aa);
        frame(2) = aa;
        bb = Basics::min(upperBound, bb);
        frame(3) = bb;
        frame(4) = upperBound;
        outputInfo.forEigenRangeCpuTime = (Real)(clock()-rangeStart) / (Real)CLOCKS_PER_SEC;
    }
    else {
        // the eigen-range is given; so do not use the Lanczos algorithm to determine it, and therefore opts.numIterForEigenRange = 0
        opts.numIterForEigenRange = 0;
        // the window of desired eigenvalues is given as [frame(2), frame(3)]
        // the range of all eigenvalues is given as [frame(1), frame(4)]
        aa = frame(2);
        bb = frame(3);
        outputInfo.forEigenRangeCpuTime = 0.0;
    }

    // set the base filter and the polynomial filter
    outputInfo.intervals = PolynomialFilterInterface::setFilter(S, frame, polyDegree, baseDegree, opts.intervalOpts);
    if (PolynomialFilterInterface::filterInfo.filterOK == 0) {
        *opts.err << "FilteredLanczosEigenSolver(Vector&, Matrix&, const SparseMatrix&, mkIndex, mkIndex, FilteredLanczosOptions&): cannot get the filter specified; please adjust your filter parameters (e.g. increasing the polynomial degree)!" << endl;
        Basics::quit(1);
    }

    // compute the eigenvectors
    opts.wantEigVec = true;
    opts.eigHighCut = PolynomialFilterInterface::filterInfo.yLimit;
    Vector dummyEigVal;
    char eigPart[] = "LA";
    LanczosInfo &shallowInfo = outputInfo;
    shallowInfo = LanczosEigenSolver(dummyEigVal, eigVec, PolynomialFilterInterface::filteredSparseMatrixPolynomialVectorProduct, S.Nrows(), opts.neigWanted, eigPart, opts);

    // compute the eigenvalues
    if (!opts.wantEigVal) {
        // if only the eigenvectors / invariant subspace are wanted, we do not have to compute the eigenvalues
        return outputInfo;  // return immediately; skip the rest which computes the eigenvalues
    }

    Matrix tmp = S*eigVec;
    mkIndex neig = dummyEigVal.Length();
    mkIndex *idx = new mkIndex[neig];
    Real *ev = new Real[neig];
    mkIndex *id = idx;
    Real *vv = ev;
    for (mkIndex i=1; i<=neig; i++) {
        *id++ = i;
        *vv++ = Vector::innerProduct(eigVec.column(i), tmp.column(i));
    }

    // reject the eigenvalues not in [aa,bb]
    mkIndex neigGood = neig;
    for (mkIndex i=0; i<neigGood; i++) {
        if ((frame(1)<aa && ev[i]<aa) || (frame(4)>bb && ev[i]>bb)) {
            // ev[i] is not a wanted eigenvalue, so exclude this one
            while (--neigGood > i) {
                if (ev[neigGood]>=aa && ev[neigGood]<=bb)
                    break;
            }
            if (neigGood > i) {
                // ev[neigGood] is a wanted eigenvalue
                // replace ev[i] by ev[neigGood], and idx[i] by idx[neigGood]
                ev[i] = ev[neigGood];
                idx[i] = idx[neigGood];
            }
        }
    }

    // sort eigenvalues ev[0],ev[1],...,ev[neigGood-1] associated with indices idx[0],idx[1],...,idx[neigGood-1]
    // by selection sort
    id = idx;
    vv = ev;
    for (mkIndex i=0; i<neigGood; i++) {
        mkIndex ii = i;
        Real val = *vv;
        Real *ww = vv;
        for (mkIndex j=i+1; j<neigGood; j++) {
            if (*(++ww) < val) {
                val = *ww;
                ii = j;
            }
        }
        if (ii != i) {
            // swap ev[i] (as *vv) and ev[ii] (as val), and idx[i] (as *id) and idx[ii]
            ev[ii] = *vv;
            *vv = val;
            mkIndex tmp = *id;
            *id = idx[ii];
            idx[ii] = tmp;
        }
        id ++;
        vv ++;
    }
    // now ev[0],ev[1],...,ev[neigGood-1] are sorted
    // the associated indices are idx[0],idx[1],...,idx[neigGood-1]

    // sort the eigenvalues
    eigVal.set(neigGood, ev);
    eigVec = eigVec.columns(neigGood, idx);

    // free local memory
    delete [] idx;

    return outputInfo;
}
