//============================================================================
// Routines for computing extreme eigenvalues of a symmetric matrix by a
// standard Lanczos procedure.
//
// Reference:
// "A Filtered Lanczos Procedure for Extreme and Interior Eigenvalue Problems"
// H.-r. Fang and Y. Saad, University of Minnesota Technical Report, 2011.
//============================================================================

#include <stdlib.h>  // for exit, atof, atoi, atol, srand, rand, etc.
#include <string.h>  // for memcpy, strcmp, strncmp, etc.
#include <time.h>    // for time_t, time, clock_t, clock, CLOCKS_PER_SEC, etc.
#include <math.h>    // for sqrt, pow, etc.
#include <iostream>  // for cout, cerr, endl, etc. (under namespace std)

#include "matkit.h"

#include "symtrieig.h"
#include "laneig.h"

#ifdef USE_MEX
    #include "mex.h"
#endif

using std::endl;

#define USE_CSRMV
// A*v corresponds to CSCMV and v*A corresponds to CSRMV
// when A is (sparse) symmetric, then A*v equals v*A in exact arithmetic
// however CSRMV usually slightly faster than CSCMV
// "#define USE_CSRMV" invokes CSRMV

////////////////////////////////////////////////////////////////////////////

class LanczosParameters {
private:
    static const SymmetricMatrix *Aptr;
    static const SparseMatrix *Sptr;
    static SymmetricMatrix A;
    static SparseMatrix S;

public:
    static void setSymmetricMatrix(const SymmetricMatrix &A0, bool shadowCopy=true) {
        if (shadowCopy)
            Aptr = &A0;
        else {
            A = A0;
            Aptr = &A;
        }
    }
    static void setSparseMatrix(const SparseMatrix &S0, bool shadowCopy=true) {
        if (shadowCopy) {
            Sptr = &S0;
        }
        else {
            S = S0;
            Sptr = &S;
        }
    }

    static Vector SymmetricMatrixVectorMultiplication(const Vector &v) {
        return (*Aptr)*v;
    }
    static Vector SparseMatrixVectorMultiplication(const Vector &v) {
        #ifdef USE_CSRMV
            return v*(*Sptr);
        #else
            return (*Sptr)*v;
        #endif
    }

    // routines for partial reorthogonalization
    static void updateOmega(mkIndex j, Real *w, const Real *w1, const Real *w2, const Real *alp, const Real *bet, Real psi, Real vartheta);
    static Real normBound(Real normb, mkIndex k, const Real *alp, const Real *bet);
};

const SymmetricMatrix *LanczosParameters::Aptr = NULL;
const SparseMatrix *LanczosParameters::Sptr = NULL;
SymmetricMatrix LanczosParameters::A;
SparseMatrix LanczosParameters::S;



////////////////////////////////////////////////////////////////////////////

// for extreme eigenvalues of a (dense) symmetric matrix
// the same as the most general form LanczosEigenSolver(), but with the NextKrylovVector() defined as A*v
LanczosInfo LanczosEigenSolver(Vector &eigVal, Matrix &eigVec, const SymmetricMatrix &A, mkIndex neigWanted, const char eigPart[], LanczosOptions &opts) {
    LanczosParameters::setSymmetricMatrix(A);
    return LanczosEigenSolver(eigVal, eigVec, LanczosParameters::SymmetricMatrixVectorMultiplication, A.Nrows(), neigWanted, eigPart, opts);
}
// the same as the most general form LanczosEigenSolver(), but with the NextKrylovVector() defined as A*v and with the default LanczosOptions
LanczosInfo LanczosEigenSolver(Vector &eigVal, Matrix &eigVec, const SymmetricMatrix &A, mkIndex neigWanted, const char eigPart[]) {
    LanczosOptions opts;
    return LanczosEigenSolver(eigVal, eigVec, A, neigWanted, eigPart, opts);
}

// for extreme eigenvalues of a sparse (symmetric) matrix
// the same as the most general form LanczosEigenSolver(), but with the NextKrylovVector() defined as A*v
LanczosInfo LanczosEigenSolver(Vector &eigVal, Matrix &eigVec, const SparseMatrix &S, mkIndex neigWanted, const char eigPart[], LanczosOptions &opts) {
    LanczosParameters::setSparseMatrix(S);
    return LanczosEigenSolver(eigVal, eigVec, LanczosParameters::SparseMatrixVectorMultiplication, S.Nrows(), neigWanted, eigPart, opts);
}
// the same as the most general form LanczosEigenSolver(), but with the NextKrylovVector() defined as A*v and with the default LanczosOptions
LanczosInfo LanczosEigenSolver(Vector &eigVal, Matrix &eigVec, const SparseMatrix &S, mkIndex neigWanted, const char eigPart[]) {
    LanczosOptions opts;
    return LanczosEigenSolver(eigVal, eigVec, S, neigWanted, eigPart, opts);
}

// the same as the most general form LanczosEigenSolver(), but with the default LanczosOptions
LanczosInfo LanczosEigenSolver(Vector &eigVal, Matrix &eigVec, NEXT_VECTOR NextKrylovVector, mkIndex n, mkIndex neigWanted, const char eigPart[]) {
    LanczosOptions opts;
    return LanczosEigenSolver(eigVal, eigVec, NextKrylovVector, n, neigWanted, eigPart, opts);
}

// Lanczos eigensolver, the most general form
// input:
// NextKrylovVector is the function which defines a self-adjoint operator, e.g. A*v with A a symmetric matrix
// n is the length of the input/output vectors of NextKrylovVector()
// neigWanted is the number of eigenvalues sought
// eigPart is a string "SA", "LA", or "BE"
//    "LA" - largest algebraic, for largest eigenvalues
//    "SA" - smallest algebraic, for smallest eigenvalues
//    "BE" - both ends, one more from high end if nev is odd
// opts is a collection of Lanczos options
// output:
// eigVal is a vector of length neigWanted containing the computed eigenvalues
// eigVec is a matrix whose columns are eigenvectors, if opts.wantEigVec==true
// return:
// this routine returns an instance of LanczosInfo which gives the information of the Lanczos procedure
LanczosInfo LanczosEigenSolver(Vector &eigVal, Matrix &eigVec, NEXT_VECTOR NextKrylovVector, mkIndex n, mkIndex neigWanted, const char eigPart[], LanczosOptions &opts) {
    // mark the time to start
    clock_t start = clock();
    clock_t forNextKrylovVectorClockTicks = 0;
    clock_t reorthogonalizationClockTicks = 0;
    clock_t convergenceCheckClockTicks = 0;
    clock_t forEigenVectorsClockTicks = 0;

    if (neigWanted == 0u) {
        // number of eigenvalues desired is not specified
        // it means that a range should given and to compute `all' eigenvalues in the range
        // check whether a range is given
        if (!strcmp(eigPart, "LA") || !strcmp(eigPart, "la")) {
            if (opts.eigHighCut == -Basics::infinity) {
                *opts.err << "LanczosEigenSolver(Vector &, Matrix &, NEXT_VECTOR, mkIndex, mkIndex, char [], LanczosOptions&): the number of eigenvalues and the lower bound are not specified (eigen part is \"LA\")!" << endl;
                Basics::quit(1);
            }
        }
        else if (!strcmp(eigPart, "SA") || !strcmp(eigPart, "sa")) {
            if (opts.eigLowCut == Basics::infinity) {
                *opts.err << "LanczosEigenSolver(Vector &, Matrix &, NEXT_VECTOR, mkIndex, mkIndex, char [], LanczosOptions&): the number of eigenvalues and the upper bound are not specified (eigen part is \"SA\")!" << endl;
                Basics::quit(1);
            }
        }
        else if (!strcmp(eigPart, "BE") || !strcmp(eigPart, "be")) {
            if (opts.eigHighCut == -Basics::infinity || opts.eigLowCut == Basics::infinity) {
                *opts.err << "LanczosEigenSolver(Vector &, Matrix &, NEXT_VECTOR, mkIndex, mkIndex, char [], LanczosOptions&): the number of eigenvalues and ";
                if (opts.eigHighCut == -Basics::infinity && opts.eigLowCut == Basics::infinity)
                    *opts.err << "the bounds";
                else if (opts.eigHighCut == -Basics::infinity)
                    *opts.err << "lower bound";
                else  // opts.eigLowCut == Basics::infinity
                    *opts.err << "upper bound";
                *opts.err << " are not specified (eigen part is \"BE\")!" << endl;
                Basics::quit(1);
            }
        }
        else {
            *opts.err << "LanczosEigenSolver(Vector &, Matrix &, NEXT_VECTOR, mkIndex, mkIndex, const char eigPart[], LanczosOption): invalid eigPart[]==\"" << eigPart << "\"!" << endl;
            Basics::quit(1);
        }
        neigWanted = n;
        if (opts.minIter == 0u)
            opts.minIter = 20;
        opts.maxIter = n;
    }

    // neigWanted, the number of total eigenvalues, is at most n
    // however, when no reorthogonalization scheme is applied, eigenvalues may re-appear
    if (opts.reorth)
        neigWanted = Basics::min(neigWanted, n);

    // preprocessing minIter, maxIter, and stride
    if (opts.minIter == 0u) {
        // set the default value
        if (!strcmp(eigPart, "BE") || !strcmp(eigPart, "be"))
            opts.minIter = opts.defaultMinIterFactor*((neigWanted+1u)/2u);
        else
            opts.minIter = opts.defaultMinIterFactor*neigWanted;
        if (opts.maxIter)
            opts.minIter = Basics::min(opts.maxIter, opts.minIter);
    }
    if (opts.maxIter == 0u) {
        // set the default value
        if (!strcmp(eigPart, "BE") || !strcmp(eigPart, "be"))
            opts.maxIter = 500 + opts.defaultMaxIterFactor*((neigWanted+1u)/2u);
        else
            opts.maxIter = 500 + opts.defaultMaxIterFactor*neigWanted;
    }
    if (opts.reorth)  {
        // if a reorthogonalization strategy is incorporated, at most n iterations are allowed
        opts.minIter = Basics::min(opts.minIter, n);
        opts.maxIter = Basics::min(opts.maxIter, n);
    }
    opts.maxIter = Basics::max(opts.maxIter, opts.minIter);
    if (opts.stride == 0u)
        opts.stride = 10;  // the default value

    // opts.memoryExpansionFactor determines how much the memory should be expanded
    // when the allocated memory is not sufficient to store alpha's and beta's
    // (the elements of the symmetric tridiagonal matrix)
    if (opts.memoryExpansionFactor <= 1.0)
        opts.memoryExpansionFactor = 1.2;
    Matrix::setMemoryExpansionFactor(opts.memoryExpansionFactor);

    // tolerance of eigenvalues and variables for convergence check
    if (opts.tol < 0.0)
        opts.tol = sqrt(Basics::machineEpsilon);
    if (opts.disp >= 2)
        *opts.out << "minIter=" << opts.minIter << ", maxIter=" << opts.maxIter << ", stride=" << opts.stride << ", extraIter=" << opts.extraIter << ", tol=" << opts.tol << endl;
    Vector ritzVal, ritzVal0;  // eigenvalues of T_k and T_{k-1}, respectively
    Real eigLowSum, eigLowSum0, absEigLowSum0;
    Real eigHighSum, eigHighSum0, absEigHighSum0;
    mkIndex neigLow = 0, neigLowFound = 0, neigLowConverged = 0;     // numbers of eigenvalues desired / found / converged in the low end of the spectrum, respectively
    mkIndex neigHigh = 0, neigHighFound = 0, neigHighConverged = 0;  // numbers of eigenvalues desired / found / converged in the high end of the spectrum, respectively
    // set the "objective"; the objective may not be met due to eigHighCut and eigLowCut
    if (!strcmp(eigPart, "LA") || !strcmp(eigPart, "la")) {
        neigHigh = neigWanted;
        neigLow = 0;
        if (opts.eigSort == -1)
            opts.eigSort = 2;  // decreasing
    }
    else if (!strcmp(eigPart, "SA") || !strcmp(eigPart, "sa")) {
        neigLow = neigWanted;
        neigHigh = 0;
        if (opts.eigSort == -1)
            opts.eigSort = 0;  // increasing
    }
    else if (!strcmp(eigPart, "BE") || !strcmp(eigPart, "be")) {
        neigLow = neigWanted / 2;         // number of desired eigenvalues at the lower end
        neigHigh = neigWanted - neigLow;  // number of desired eigenvalues at the higher end (one more than neigLow if neigWanted is odd)
        if (opts.eigSort == -1)
            opts.eigSort = 0;  // increasing both ends
    }
    else {
        *opts.err << "LanczosEigenSolver(Vector &, Matrix &, NEXT_VECTOR, mkIndex, mkIndex, const char eigPart[], LanczosOption): invalid eigPart[]==\"" << eigPart << "\"!" << endl;
        Basics::quit(1);
    }
    // output information
    LanczosInfo outputInfo;

    // initialization
    Vector v0(n), v1, v2;
    if (opts.v0.Length() == n)
        v1 = opts.v0;
    else {
        if (opts.disp && opts.v0.Length()>0u)
            *opts.out << "laneig warning: the initial Lanczos vector is of size " << opts.v0.Length() << ", not of size " << n << "; will use a random vector instead!" << endl;
        v1 = Vector::random(n, -1.0, 1.0);  // another option is v1 = Vector::Ones(n); however it may produce poor results in some cases
    }
    v1 /= v1.norm2();
    Vector *pv0 = &v0;
    Vector *pv1 = &v1;
    Vector *pv2 = &v2;
    Matrix V(n, 0, opts.minIter);
    mkIndex memorySize = opts.minIter;
    Real *alp = new Real[memorySize];     // diagonal of the Lanczos matrix T
    Real *bet = new Real[memorySize];     // off-diagonal of the Lanczos matrix T
    Real aa, bb_old, bb = 0.0;            // aa and bb are the alpha and beta at the current Lanczos iteration; bb_old is the previous beta
                                          // bb is initialized as 0.0 for the first iteration
    Real betaSum = 0.0;                   // the sum of all beta's

    // variables for reorthogonalization count
    mkIndex reorthCountThisIteration = 0;
    // variables for partial reorthogonalization
    bool forceReorth = false;
    // eps1 is a bound, in practice, of inner product of two numerically orthogonal vectors of length n
    // note that machineEpsilon/2.0 is the unit roundoff
    const Real eps1 = sqrt((Real)n)*Basics::machineEpsilon/2.0;
    // two vectors are called semi-orthogonal to each other if their inner product is <= delta
    const Real delta = sqrt(Basics::machineEpsilon);
    // when the estimated bound on inner product of two Lanczos vectors < eta, reorthogonalization is performed
    const Real eta = Basics::min(eps1, (Real)pow(Basics::machineEpsilon, 3.0/4.0));
    // normb is a bound on ||T_k||_2, where T_k is the k-by-k symmetric tridiagonal Lanczos matrix
    Real normb = 0.0;
    Real *w0 = NULL, *w1 = NULL, *w2 = NULL;
    bool *flag = NULL;
    if (opts.reorth == 1) {
        // partial reorthogonalization
        w0 = new Real[opts.maxIter];
        w1 = new Real[opts.maxIter];
        w2 = new Real[opts.maxIter];
        w1[0] = 1.0;  // in the second iteration, it==1, pointer w2 points to the initial w1, and
                      // initial w1[0] is v'*v, where v is the first Lanczos vector
        if (opts.partialReorthStrategy <= 2)
            flag = new bool[opts.maxIter];
        outputInfo.maxOrthLevelRatio = 0.0;
        outputInfo.minOrthLevelRatio = Basics::infinity;
    }

    // the big Lanczos loop
    mkIndex iter = 0;
    outputInfo.allEigenvaluesCheckedConverged = false;
    outputInfo.reorthVectorCount = 0;
    outputInfo.reorthIterCount = 0;
    outputInfo.doubleReorthIterCount = 0;
    outputInfo.localReorthIterCount = 0;
    mkIndex minIter = opts.minIter;
    while (true) {
        if (iter == memorySize) {
            if (opts.disp >= 2) {
                *opts.out << "at iteration #" << iter+1u << ", memory expanded for storing Lanczos vectors and alphas and betas" << endl;
            }
            // allocated memory is insufficient, so expand
            memorySize *= opts.memoryExpansionFactor;
            if (memorySize <= iter)
                memorySize = iter+1u;
            Real *tmp = new Real[memorySize];
            memcpy(tmp, alp, iter*sizeof(Real));
            delete [] alp;
            alp = tmp;
            tmp = new Real[memorySize];
            memcpy(tmp, bet, iter*sizeof(Real));
            delete [] bet;
            bet = tmp;
        }
        V.appendColumn(pv1->Store());   // append the new Lanczos vector v1 to matrix V
                                        // the memory expansion is also applied here to store Lanczos vectors, if iter == memorySize

        clock_t forNextKrylovVectorStart = clock();
        *pv2 = NextKrylovVector(*pv1);  // the next Lanczos vector before orthogonalization
        forNextKrylovVectorClockTicks += (clock()-forNextKrylovVectorStart);
        (*pv2) -= bb*(*pv0);            // orthogonalization against v0, the vector before the previous vector v1
        // the above code essentially computes v2 = A*v1 - bb*v0, but we replace A*v1 by NextKrylovVector(v1) for general cases

        aa = Vector::innerProduct(*pv1, *pv2);  // aa = <v1,v2>  (the new alpha)
        alp[iter] = aa;                 // add the diagonal element to the Lanczos matrix T

        // a special case: max # iterations reached, quit anyway (do not have to check convergence)
        if (iter+1u == opts.maxIter) {
            if (opts.wantEigVec) {
                (*pv2) -= aa*(*pv1);    // reorthogonalization  (for the last beta)
                bet[iter] = pv2->norm2();
                // the last beta, coupled with the bottom elements of Ritz vectors (i.e. eigenvectors of T), can be used for computing eigenvalue errors
                // see Parlett's book "The Symmetric Eigenvalue Problem", page 290
            }
            iter ++;
            break;
        }

        // check convergence
        if ((opts.reorth==0 || iter+1u<n) && iter+1u>=minIter && (iter+1u-minIter)%opts.stride==0u) {
            clock_t checkStart = clock();
            // the last beta, coupled with the bottom elements of Ritz vectors (i.e. eigenvectors of T), can be used for computing eigenvalue errors
            // see Parlett's book "The Symmetric Eigenvalue Problem", page 290
            int info;
            if (opts.stride>1 || iter+1u==minIter) {
                info = SymmetricTridiagoanlEigenSolver(ritzVal0, iter, alp, bet);
                if (opts.disp)
                    reportTroubleIfAny(*opts.out, info, iter);
            }
            info = SymmetricTridiagoanlEigenSolver(ritzVal, iter+1u, alp, bet);
            if (opts.disp)
                reportTroubleIfAny(*opts.out, info, iter+1u);
            // check the convergence of eigenvalues at the low end
            mkIndex nev = neigLow;  // number of smallest eigenvalues wanted
            const Real *rv0 = ritzVal0.Store();
            const Real *rv = ritzVal.Store();
            eigLowSum = 0.0;
            eigLowSum0 = 0.0;
            absEigLowSum0 = 0.0;
            while (nev && (*rv) <= opts.eigLowCut) {
                absEigLowSum0 += Basics::abs(*rv0);
                eigLowSum0 += (*rv0++);
                eigLowSum += (*rv++);
                nev --;
            }
            neigLowFound = neigLow - nev;
            Real eigLowErr = Basics::abs((eigLowSum-eigLowSum0) / absEigLowSum0);
            if (neigLow && opts.disp >= 2)
                *opts.out << "at iteration #" << iter+1u << ", estimated eigen error (low end) is " << eigLowErr << " (" << neigLow-nev << " eigenvalues)" << endl;
            if (absEigLowSum0 == 0.0 && eigLowSum == 0.0)
                eigLowErr = 0.0;  // a special case, considered as no error
            // check the convergence of eigenvalues at the high end
            mkIndex nev2 = neigHigh;
            rv0 = ritzVal0.Store() + iter - 1;
            rv = ritzVal.Store() + iter;
            eigHighSum = 0.0;
            eigHighSum0 = 0.0;
            absEigHighSum0 = 0.0;
            while (nev2 && (*rv) >= opts.eigHighCut) {
                absEigHighSum0 += Basics::abs(*rv0);
                eigHighSum0 += (*rv0--);
                eigHighSum += (*rv--);
                nev2 --;
            }
            neigHighFound = neigHigh - nev2;
            Real eigHighErr = Basics::abs((eigHighSum-eigHighSum0) / absEigHighSum0);
            if (neigHigh && opts.disp >= 2) {
                *opts.out << "at iteration #" << iter+1u << ", estimated eigen error (high end) is " << eigHighErr << " (" << neigHigh-nev2 << " eigenvalues)" << endl;
            }
            if (absEigHighSum0 == 0.0 && eigHighSum == 0.0)
                eigHighErr = 0.0;  // a special case, considered as no error
            if (eigLowErr <= opts.tol && eigHighErr <= opts.tol) {
                if (opts.extraIter == 0u) {
                    neigLowConverged = neigLowFound;
                    neigHighConverged = neigHighFound;
                }
                if (neigLowConverged < neigLowFound || neigHighConverged < neigHighFound) {
                    neigLowConverged = neigLowFound;
                    neigHighConverged = neigHighFound;
                    minIter = Basics::min(iter+1+opts.extraIter, opts.maxIter);
                    if (opts.disp >= 2) {
                        *opts.out << "at iteration #" << iter+1u << ", eigenvalues are deemed converged, now running extra iterations" << endl;
                    }
                }
                else {
                    outputInfo.allEigenvaluesCheckedConverged = true;
                    convergenceCheckClockTicks += (clock()-checkStart);
                    if (opts.wantEigVec) {
                        (*pv2) -= aa*(*pv1);
                        bb = pv2->norm2();
                        bet[iter] = bb;
                        // the last beta, coupled with the bottom elements of Ritz vectors (i.e. eigenvectors of T), can be used for computing eigenvalue errors
                        // see Parlett's book "The Symmetric Eigenvalue Problem", page 290
                    }
                    iter ++;
                    break;  // all wanted eigenvalues converged; quit the big Lanczos loop
                }
            }
            if (opts.stride == 1)
                ritzVal0 = ritzVal;
            convergenceCheckClockTicks += (clock()-checkStart);
        }  // end of check convergence

        (*pv2) -= aa*(*pv1);            // v2 = v2 - aa*v1, reorthogonalization against the previous Lanczos vector v1, subject to reorthogonalization and normalization
        bb = pv2->norm2();              // the new beta
        bet[iter] = bb;                 // add the new off-diagonal element to the symmetric tridiagonal matrix T
        betaSum += bb;

        // reorthogonalization
        clock_t reorthogonalizationStart = clock();
        if (opts.reorth == 2) {  // full reorthogonalization
            // to reorthogonalize against all previous iter+1 Lanczos vectors
            outputInfo.reorthIterCount ++;
            outputInfo.reorthVectorCount += (iter+1u);
            pv2->modifiedGramSchmidt(V);
            // the following two lines are the classical Gram-Schmidt, whereas the above line is the modified Gram-Schmidt
            // Vector c = (*pv2)*V;     // c = (v2'*V)'
            // (*pv2) -= V*c;           // v2 = v2 - V*c
            bb_old = bb;
            bb = pv2->norm2();
            if (opts.doubleReorthGamma >= 1.0 || bb < opts.doubleReorthGamma*bb_old) {
                outputInfo.doubleReorthIterCount ++;
                if (opts.disp >= 2 && opts.doubleReorthGamma < 1.0) {
                    *opts.out << "at iteration #" << iter+1u << ", reorthogonalization is doubled" << endl;
                }
                outputInfo.reorthVectorCount += (iter+1u);
                pv2->modifiedGramSchmidt(V);
                bb = pv2->norm2();
            }
            bet[iter] = bb;  // update the latest beta
        }
        else if (opts.reorth == 1) {  // partial reorthogonalization
            // extended local reorthogonalization
            if (opts.localReorthGamma >= 1.0 || (iter >= 1u && bb < opts.localReorthGamma*bet[iter-1u])) {
                if (opts.disp >= 2 && opts.localReorthGamma < 1.0) {
                    *opts.out << "at iteration #" << iter+1u << ", extended local reorthogonalization is performed" << endl;
                }
                // reorthogonalize against the Lanczos vector two iterations ago
                outputInfo.localReorthIterCount ++;
                if (iter > 0u) {
                    outputInfo.reorthVectorCount ++;
                    bb = Vector::innerProduct(*pv2, *pv0);  // orthogonalize again against v0, the vector before previous Lanczos vector
                    (*pv2) -= bb*(*pv0);
                    if (bet[iter-1u] != 0.0)  // reduce risk of cancellation
                        bet[iter-1u] += bb;
                }
                // reorthogonalize against the previous Lanczos vector
                outputInfo.reorthVectorCount ++;
                aa = Vector::innerProduct(*pv2, *pv1);
                (*pv2) -= aa*(*pv1);
                alp[iter] += aa;
                bb = pv2->norm2();
                bet[iter] = bb;  // update the latest beta, which will be used in the routine updateOmega for partial reorthogonalization
            }
            // update estimated level of orthogonality
            normb = LanczosParameters::normBound(normb, iter+1u, alp, bet);
            LanczosParameters::updateOmega(iter, w0, w1, w2, alp, bet, eps1, eps1*normb);  // update omega, an estimated upper bound of v_i'*v_j
            Real max_w = 0.0;
            bool *flagx = NULL;
            if (forceReorth) {
                // reorthogonalize for semi-orthogonality
                mkIndex ii = iter;
                Real *ww = w0;
                while (ii--)
                    *(ww++) = eps1;
                reorthCountThisIteration = iter+1;  // to reorthogonalize against all previous iter+1 Lanczos vectors
                flagx = NULL;                       // means all previous Lanczos vectors are flagged to reorthogonalize against to
            }
            else {
                const Real *ww = w0;
                max_w = Basics::abs(*ww++);
                mkIndex ii = iter;
                while (ii--)
                    max_w = Basics::max(max_w, Basics::abs(*ww++));
                // check reorth if requested, the cost of checking is high, about half of the full reorthogonalization
                if (opts.checkReorth) {
                    Vector c = ((*pv2)*V)/bb;  // equivalently, c = (v2'*V)', with v2 being normalized;
                                               // c(i) is the inner product of the current (last) Lanczos vector and the i-th Lanczos vector
                    Real max_w0 = 0.0;
                    ii = c.Length();
                    ww = c.Store();
                    while (ii--)
                        max_w0 = Basics::max(max_w0, Basics::abs(*ww++));
                    if (max_w0 > max_w) {
                        if (max_w <= delta && max_w0 > delta)
                            *opts.out << "laneig error: at iteration #" << iter+1u << ", losing semi-orthogonality (tol=" << delta << ") because estimated max_w (" << max_w << ") is less than true max_w (" << max_w0 << ")!" << endl;
                        else
                            *opts.out << "laneig warning: at iteration #" << iter+1u << ", estimated max_w (" << max_w << ") is less than true max_w (" << max_w0 << ")!" << endl;
                    }
                    outputInfo.minOrthLevelRatio = Basics::min(outputInfo.minOrthLevelRatio, max_w/max_w0);
                    outputInfo.maxOrthLevelRatio = Basics::max(outputInfo.maxOrthLevelRatio, max_w/max_w0);
                }
                if (max_w > delta) {
                    // reorthogonalize for semi-orthogonality
                    if (opts.partialReorthStrategy == 0) {
                        // cheapest
                        reorthCountThisIteration = 0;
                        bool *flg = flag;
                        mkIndex ii = iter+1u;
                        // reset flags
                        while (ii--)
                            *flg++ = false;
                        Real *ww = w0;
                        for (ii=0; ii<=iter; ii++) {
                            if (Basics::abs(*ww) > delta) {
                                Real *w2 = ww-1;
                                flg = flag+ii;
                                *flg = true;
                                // look backward for abs(w0[k]) > eta
                                for (mkIndex jj=ii; jj>0; jj--) {
                                    flg --;
                                    if (Basics::abs(*w2--) > eta)
                                        *flg = true;
                                }
                                w2 = ww+1;
                                flg = flag+ii;
                                // look forward for abs(w0[k]) > eta
                                for (mkIndex jj=ii; jj<iter; jj++) {
                                    flg ++;
                                    if (Basics::abs(*w2++) > eta)
                                        *flg = true;
                                }
                            }
                            ww++;
                        }
                        ii = iter+1;
                        ww = w0;
                        flg = flag;
                        while (ii--) {
                            if (*flg++) {
                                *ww = eps1;
                                reorthCountThisIteration ++;
                            }
                            ww++;
                        }
                        flagx = flag;
                    }
                    else if (opts.partialReorthStrategy == 1) {
                        // second to the cheapest
                        reorthCountThisIteration = 0;
                        bool *flg = flag;
                        Real *ww = w0;
                        mkIndex ii = iter+1u;
                        while (ii--) {
                            if (Basics::abs(*ww) <= eta) {
                                *flg++ = false;  // not to reorthogonalize against to
                                ww++;
                            }
                            else {
                                *flg++ = true;   // to reorthogonalize against to
                                *ww++ = eps1;
                                reorthCountThisIteration ++;
                            }
                        }
                        flagx = flag;
                    }
                    else if (opts.partialReorthStrategy == 2) {
                        // more expensive
                        bool *flg = flag;
                        Real *ww = w0;
                        mkIndex ii = 0;
                        for (; ii<=iter; ii++) {
                            if (Basics::abs(*ww++) <= eta)
                                *flg++ = false;  // not to reorthogonalize against to
                            else
                                break;
                        }
                        mkIndex jj = iter;
                        flg = flag + iter;
                        ww = w0 + iter;
                        for (; jj>ii; jj--) {
                            if (Basics::abs(*ww--) <= eta)
                                *flg-- = false;  // not to reorthogonalize against to
                            else
                                break;
                        }
                        ww = w0 + ii;
                        flg = flag + ii;
                        mkIndex steps = jj+1u-ii;
                        reorthCountThisIteration = steps;
                        while (steps--) {
                            *flg++ = true;       // to reorthogonalize against to
                            *ww++ = eps1;
                        }
                        flagx = flag;
                    }
                    else {  // opts.partialReorthStrategy >= 3
                        // reorthogonalize against all previous vectors, most expensive
                        mkIndex ii = iter+1u;
                        Real *ww = w0;
                        while (ii--)
                            *ww++ = eps1;
                        // the w0[iter] was initialized as eps1, so it is not necessary set it as eps1 here
                        reorthCountThisIteration = iter+1;  // to reorthogonalize against all previous iter+1 Lanczos vectors
                        flagx = NULL;                       // means all previous Lanczos vectors are flagged to reorthogonalize against to
                    }
                }
            }
            if (forceReorth || max_w>delta) {
                // perform partial reorthogonalization
                if (opts.disp >= 2) {
                    if (forceReorth)
                        *opts.out << "at iteration #" << iter+1u << ", reorthogonalization is enforced" << endl;
                    else
                        *opts.out << "at iteration #" << iter+1u << ", reorthogonalization is performed" << endl;
                }
                outputInfo.reorthIterCount ++;
                outputInfo.reorthVectorCount += reorthCountThisIteration;
                pv2->modifiedGramSchmidt(V, flagx);
                bb_old = bb;
                bb = pv2->norm2();
                if (bb < opts.doubleReorthGamma*bb_old) {
                    outputInfo.doubleReorthIterCount ++;
                    if (opts.disp >= 2) {
                        *opts.out << "at iteration #" << iter+1u << ", reorthogonalization is doubled" << endl;
                    }
                    outputInfo.reorthVectorCount += reorthCountThisIteration;
                    pv2->modifiedGramSchmidt(V, flagx);
                    bb = pv2->norm2();
                }
                bet[iter] = bb;  // update latest beta
                forceReorth = !forceReorth;
            }
            Real *tmp = w2;
            w2 = w1;
            w1 = w0;
            w0 = tmp;
        }
        if (bb < betaSum*Basics::machineEpsilon || bb == 0.0) {
            // special case, the new beta (bb) is deemed 0
            (*pv2) = Vector::random(n, -1.0, 1.0);
            pv2->modifiedGramSchmidt(V);
            bb = pv2->norm2();
            bet[iter] = 0.0;
            if (opts.reorth == 1) {
                // partial reorthogonalization
                forceReorth = true;
            }
        }
        reorthogonalizationClockTicks += (clock()-reorthogonalizationStart);
        *pv2 /= bb;              // v2 = v2 / bet[iter], except for the special case

        // update pointers pv1 <- pv2, pv2 <- pv0, pv0 <- pv1
        Vector *tmp = pv1;
        pv1 = pv2;               // so now *pv1 is the current Lanczos vector and
        pv2 = pv0;               // *pv2 is for the next Lanczos vector and
        pv0 = tmp;               // *pv0 is the previous Lanczos vector
        iter ++;
        // the following lines print out some information of the level of orthogonality (expensive)
        // *opts.out << "||V'*V-I||=" << (V.t()*V-SparseMatrix::Eye(iter+1)).normFrobenius() << endl;
        // SparseMatrix T = SparseMatrix::tridiagonalMatrix(iter+1,bet,alp,bet);
        // *opts.out << "||V'*A*V-T||=" << (V.t()*(A*V)-T).normFrobenius() << endl;
        // *opts.out << "||V*T*V'-A||=" << (V*(T*V.t())-A).normFrobenius() << endl;
    }   // end of the big Lanczos loop

    // free some memory that was used for partial reorthogonalization
    if (opts.reorth == 1) {
        delete [] w0;
        delete [] w1;
        delete [] w2;
        if (opts.partialReorthStrategy <= 2)
            delete [] flag;
    }

    // number of Lanczos iterations used, and number of vectors against which the reorthogonalization is performed
    outputInfo.numIter = iter;
    if (iter == 1u)
        outputInfo.reorthVectorRate = 0.0;
    else
        outputInfo.reorthVectorRate = (((Real)outputInfo.reorthVectorCount/(Real)iter)/(Real)(iter-1u))*2.0;
    if (opts.disp && opts.reorth==1 && opts.checkReorth) {
        *opts.out << "max orthogonality level ratio: " << outputInfo.maxOrthLevelRatio << endl;
        *opts.out << "min orthogonality level ratio: " << outputInfo.minOrthLevelRatio << endl;
    }
    if (opts.disp >= 2) {
        *opts.out << "reorthogonalization vector count: " << outputInfo.reorthVectorCount << " (" << 100.0*outputInfo.reorthVectorRate << "%%)" << endl;
    }

    Matrix S;
    int info = 0;
    if (opts.wantEigVec) {
        clock_t eigvecStart = clock();
        info = SymmetricTridiagoanlEigenSolver(ritzVal, S, iter, alp, bet);
        forEigenVectorsClockTicks += clock()-eigvecStart;
    }
    else if (ritzVal.Length() != iter) {
        // maximum number of iterations is met
        info = SymmetricTridiagoanlEigenSolver(ritzVal, iter, alp, bet);
    }
    if (opts.disp)
        reportTroubleIfAny(*opts.out, info);

    if (!outputInfo.allEigenvaluesCheckedConverged) {
        if (opts.disp && (opts.reorth==0 || opts.maxIter<n) && opts.minIter!=opts.maxIter) {
            *opts.out << "laneig warning: the maximum number of Lanczos iterations " << iter << " is met, but the desired eigenvalues may not yet be all converged!"<< endl;
        }
        const Real *rv = ritzVal.Store();
        mkIndex nev = neigLow;           // number of smallest eigenvalues wanted
        while (nev && (*rv++) <= opts.eigLowCut)
            nev --;
        neigLowFound = neigLow - nev;    // number of smallest eigenvalues found
        rv = ritzVal.Store() + iter - 1;
        nev = neigHigh;                  // number of largest eigenvalues wanted
        while (nev && (*rv--) >= opts.eigHighCut)
            nev --;
        neigHighFound = neigHigh - nev;  // number of largest eigenvalues found
    }

    // neigLowFound is the number of smallest eigenvalues found; neigHighFound is the number of largest eigenvalues found
    mkIndex a, b, c, d;
    if (opts.eigSort == 1) {
        a = neigLowFound;
        b = 1;
        c = iter-neigHighFound+1;
        d = iter;
    }
    else if (opts.eigSort == 2) {
        a = 1;
        b = neigLowFound;
        c = iter;
        d = iter-neigHighFound+1;
    }
    else if (opts.eigSort == 3) {
        a = neigLowFound;
        b = 1;
        c = iter;
        d = iter-neigHighFound+1;
    }
    else if (opts.eigSort == 4) {
        a = iter-neigHighFound+1;
        b = iter;
        c = 1;
        d = neigLowFound;
    }
    else if (opts.eigSort == 5) {
        a = iter-neigHighFound+1;
        b = iter;
        c = neigLowFound;
        d = 1;
    }
    else if (opts.eigSort == 6) {
        a = iter;
        b = iter-neigHighFound+1;
        c = 1;
        d = neigLowFound;
    }
    else if (opts.eigSort == 7) {
        a = iter;
        b = iter-neigHighFound+1;
        c = neigLowFound;
        d = 1;
    }
    else {
        // opts.eigSort == 0 (increasing) || opts.eigSort <= -2 (no sort)
        // the case opts.eigSort == -1 (default) should have been handled
        a = 1;
        b = neigLowFound;
        c = iter-neigHighFound+1;
        d = iter;
    }

    // now get eigenvalues from the two ends of spectrum
    Matrix bottomElements;
    if (neigHighFound == 0u) {
        if (neigLowFound == 0u) {
            eigVal = Vector();
            if (opts.wantEigVec)
                eigVec = Matrix();
        }
        else {
            eigVal = ritzVal.subVector(a,b);
            if (opts.wantEigVec) {
                // computing eigenvectors
                clock_t eigvecStart = clock();
                eigVec = V*S.columns(a,b);
                forEigenVectorsClockTicks += clock()-eigvecStart;
                bottomElements = S.subMatrix(iter,iter,a,b);
            }
        }
    }
    else if (neigLowFound == 0u) {
        eigVal = ritzVal.subVector(c,d);
        if (opts.wantEigVec) {
            // computing eigenvectors
            clock_t eigvecStart = clock();
            eigVec = V*S.columns(c,d);
            forEigenVectorsClockTicks += clock()-eigvecStart;
            bottomElements = S.subMatrix(iter,iter,c,d);
        }
    }
    else {
        eigVal = Vector::concatenate(ritzVal.subVector(a,b), ritzVal.subVector(c,d));
        if (opts.wantEigVec) {
            // computing eigenvectors
            clock_t eigvecStart = clock();
            eigVec = V * Matrix::horizontalConcatenate(S.columns(a,b), S.columns(c,d));
            forEigenVectorsClockTicks += clock()-eigvecStart;
            bottomElements = Matrix::horizontalConcatenate(S.subMatrix(iter,iter,a,b), S.subMatrix(iter,iter,c,d));
        }
    }

    // compute eigenvectors, if desired
    if (opts.wantEigVec) {
        // the last beta, coupled with the bottom elements of Ritz vectors (i.e. eigenvectors of T), can be used for computing eigenvalue errors
        // see Parlett's book "The Symmetric Eigenvalue Problem", page 290
        Real maxErr0 = 0.0, maxErr = 0.0;
        Real bottomBeta = bet[iter-1u];
        for (mkIndex i=1; i<=neigLowFound+neigHighFound; i++) {
            Real err0 = Basics::abs(bottomBeta*bottomElements(1,i));
            // abs(bet[iter-1]*S(iter,i)) is identical to (A*x-ritzVal(i)*x).norm2() subject to rounding errors
            maxErr0 = Basics::max(maxErr0, err0);
            maxErr = Basics::max(maxErr, Basics::abs(err0/ritzVal(i)));
        }
        outputInfo.maxEigenvalueError = maxErr0;
        outputInfo.maxRelativeEigenvalueError = maxErr;
    }
    delete [] alp;
    delete [] bet;

    // CPU time
    outputInfo.LanczosCpuTime = (Real)(clock()-start) / (Real)CLOCKS_PER_SEC;
    outputInfo.reorthogonalizationCpuTime = (Real)reorthogonalizationClockTicks / (Real)CLOCKS_PER_SEC;
    outputInfo.forNextKrylovVectorCpuTime = (Real)forNextKrylovVectorClockTicks / (Real)CLOCKS_PER_SEC;
    outputInfo.convergenceCheckCpuTime = (Real)convergenceCheckClockTicks / (Real)CLOCKS_PER_SEC;
    outputInfo.forEigenVectorsCpuTime = (Real)forEigenVectorsClockTicks / (Real)CLOCKS_PER_SEC;

    // required memory
    outputInfo.memoryForLanczosInBytes = (n*V.Maxncols())*sizeof(Real);  // for storing Lanczos vectors
    outputInfo.memoryForLanczosInBytes += (5*sizeof(Real)+sizeof(bool))*memorySize;  // for alp, bet, w0, w1, w2, flag
    return outputInfo;
}


// w_{j,k} is the estimated upper bound of v_j^T*v_k to measure the satisfaction of semi-orthogonality
// assuming w_{j,*} and w_{j-1,*} are available as input w1 and w2, compute w_{j+1,*} as output w
void LanczosParameters::updateOmega(mkIndex j, Real *w, const Real *w1, const Real *w2, const Real *alp, const Real *bet, Real psi, Real vartheta) {
    // a special case j==0
    if (j == 0u) {
        w[0] = psi;
        w[1] = 1.0;
        return;
    }

    // some variables
    const Real alpj = alp[j], betj = bet[j], betj1 = bet[j-1u];
    Real val;

    // compute w[k] with k==0
    val = bet[0]*w1[1] - betj1*w2[0] + (alp[0]-alpj)*w1[0];
    if (val >= 0.0)
        val += vartheta;
    else
        val -= vartheta;
    w[0] = val/betj;

    // compute w[k] for k=1,...,j-1
    const Real *alpk = alp+1;
    const Real *betk = bet+1;
    const Real *w1k = w1+1;
    const Real *w2k = w2+1;
    Real *wk = w+1;
    mkIndex kk = j-1u;
    while (kk--) {
        // compute w_{j+1,k} = (1/beta_{j+1}) * ( beta_{k+1}*w_{j,k+1} + (alpha_k-alpha_j)*w_{j,k} + beta_k*w_{j,k-1} - beta_j*w_{j-1,k} + vartheta_{j,k} )
        // note that bet[k] is beta_{k+1} in Simon's papers (1984)
        //           input  w_{j,*}   is stored in w1
        //           input  w_{j-1,*} is stored in w2
        //           output w_{j+1,*} is stored in w
        val = (*betk)*(*(w1k+1)) - betj1*(*w2k) + ((*alpk)-alpj)*(*w1k) + *(betk-1)*(*(w1k-1));
        alpk ++;
        betk ++;
        w1k ++;
        w2k ++;
        // in other words, val = bet[k]*w1[k+1] - bet[j+1]*w2[k] + (alp[k]-alp[j])*w1[k] + bet[k-1]*w1[k-1];
        // where kk+k==j before "kk--" (or kk+k==j-1 after "kk--")
        if (val >= 0.0)
            val += vartheta;
        else
            val -= vartheta;
        *(wk++) = val/betj;
    }

    // compute w[k] for k==j and k=j+1
    *(wk++) = psi;  // w[j] = psi;
    *wk = 1.0;      // w[j+1] = 1.0;
}

// normBound provides an upper bound of ||T_k||_2, where T_k is a
// symmetric tridiagonal matrix with diagonal elements are
// alp[0],...,alp[k-1] and off-diagonal elements bet[0],...,bet[k-2].
// ||T_k||_2 is square root of ||T_k*T_k||_2 which is symmetric, so
// its 2-norm is the largest eigenvalue in magnitude which has an upper
// bound by the Gershgorin circle theorem.
// T_k*T_k is denoted by T2 in the code.
Real LanczosParameters::normBound(Real normb, mkIndex k, const Real *alp, const Real *bet) {
    Real nrm2;
    static Real bound2, bound1;
    // the rows/columns of T_k*T_k are indexed by 0,1,...,k-1
    // bound2 is the Gershgorin circle bound of row/column k-3 of T_{k-1}*T_{k-1}
    // bound1 is the Gershgorin circle bound of row/column k-2 of T_{k-1}*T_{k-1}

    if (k == 0u)
        normb = 0.0;
    else if (k == 1u)
        normb = Basics::abs(alp[0]);
    else if (k == 2u) {
        // T2 = [  alp[0]^2+bet[0]^2,     (alp[0]+alp[1])*bet[0];
        //        (alp[0]+alp[1])*bet[0],  alp[1]^2+bet[0]^2      ]
        Real radius = Basics::abs((alp[0]+alp[1])*bet[0]);
        bound1 = bet[0]*bet[0];
        bound2 = bound1;
        bound1 += (alp[1]*alp[1] + radius);
        bound2 += (alp[0]*alp[0] + radius);
        normb = sqrt(Basics::max(bound2, bound1));
    }
    else {  // k >= 3
        Real b0 = bet[k-3u], b1 = bet[k-2u];
        Real a0 = alp[k-2u], a1 = alp[k-1u];
        Real tmp0 = Basics::abs(b0*b1);
        Real tmp1 = b1*b1;
        Real tmp2 = Basics::abs((a0+a1)*b1);
        bound2 += tmp0;                       // Gershgorin circle bound for row/column k-3 of T_k*T_k
        bound1 += (tmp1 + tmp2);              // Gershgorin circle bound for row/column k-2 of T_k*T_k
        nrm2 = Basics::max(bound2, bound1);
        bound2 = bound1;
        bound1 = a1*a1 + tmp1 + tmp0 + tmp2;  // Gershgorin circle bound for row/column k-1 of T_k*T_k
        nrm2 = Basics::max(nrm2, bound1);
        normb = Basics::max(normb, (Real)sqrt(nrm2));
    }
    return normb;
}

// functionally the same as norm_bound, but not taking full advantage of incremental update
/* Real norm_bound0(Real normb, mkIndex k, const Real *alp, const Real *bet) {
    Real nrm;
    mkIndex i;

    if (k == 0u)
        normb = 0.0;
    else if (k == 1u)
        normb = Basics::abs(alp[0]);
    else if (k == 2u) {
        // T2 = [  alp[0]^2+bet[0]^2,      (alp[0]+alp[1])*bet[0];
        //         (alp[0]+alp[1])*bet[0],  alp[1]^2+bet[0]^2      ]
        Real b0 = (alp[0]+alp[1])*bet[0];
        Real a0 = bet[0]*bet[0];
        Real a1 = a0;
        a0 += alp[0]*alp[0];
        a1 += alp[1]*alp[1];
        normb = sqrt(Basics::max(a0+Basics::abs(b0), a1+Basics::abs(b0)));
    }
    else if (k == 3) {
        i = 0;
        normb = sqrt(alp[i]*alp[i] + bet[i]*bet[i] +                    // T2(i,i)
                   Basics::abs(alp[i]*bet[i] + alp[i+1]*bet[i]) +       // T2(i+1,i)
                   Basics::abs(bet[i]*bet[i+1]) );                      // T2(i+2,i)
        i = 1;
        nrm = sqrt(Basics::abs(bet[i-1]*alp[i-1] + alp[i]*bet[i-1]) +   // T2(i-1,i)
                   bet[i-1]*bet[i-1] + alp[i]*alp[i] + bet[i]*bet[i] +  // T2(i,i)
                   Basics::abs(alp[i]*bet[i] + alp[i+1]*bet[i]) );      // T2(i+1,i)
        normb = Basics::max(normb, nrm);
        i = 2;
        nrm = sqrt(Basics::abs(bet[i-1]*bet[i-2]) +                     // T2(i-2,i)
                   Basics::abs(bet[i-1]*alp[i-1] + alp[i]*bet[i-1]) +   // T2(i-1,i)
                   bet[i-1]*bet[i-1] + alp[i]*alp[i] );                 // T2(i,i)
        normb = Basics::max(normb, nrm);
    }
    else if (k == 4) {
        i = 1;
        nrm = sqrt(Basics::abs(bet[i-1]*alp[i-1] + alp[i]*bet[i-1]) +   // T2(i-1,i)
                   bet[i-1]*bet[i-1] + alp[i]*alp[i] + bet[i]*bet[i] +  // T2(i,i)
                   Basics::abs(alp[i]*bet[i] + alp[i+1]*bet[i]) +       // T2(i+1,i)
                   Basics::abs(bet[i]*bet[i+1]) );                      // T2(i+2,i)
        normb = Basics::max(normb, nrm);
        i = 2;
        nrm = sqrt(Basics::abs(bet[i-1]*bet[i-2]) +                     // T2(i-2,i)
                   Basics::abs(bet[i-1]*alp[i-1] + alp[i]*bet[i-1]) +   // T2(i-1,i)
                   bet[i-1]*bet[i-1] + alp[i]*alp[i] + bet[i]*bet[i] +  // T2(i,i)
                   Basics::abs(alp[i]*bet[i] + alp[i+1]*bet[i]) );      // T2(i+1,i)
        normb = Basics::max(normb, nrm);
        i = 3;
        nrm = sqrt(Basics::abs(bet[i-1]*bet[i-2]) +                     // T2(i-2,i)
                   Basics::abs(bet[i-1]*alp[i-1] + alp[i]*bet[i-1]) +   // T2(i-1,i)
                   bet[i-1]*bet[i-1] + alp[i]*alp[i] );                 // T2(i,i)
        normb = Basics::max(normb, nrm);
    }
    else {
        i = k-3;
        nrm = sqrt(Basics::abs(bet[i-1]*bet[i-2]) +                     // T2(i-2,i)
                   Basics::abs(bet[i-1]*alp[i-1] + alp[i]*bet[i-1]) +   // T2(i-1,i)
                   bet[i-1]*bet[i-1] + alp[i]*alp[i] + bet[i]*bet[i] +  // T2(i,i)
                   Basics::abs(alp[i]*bet[i] + alp[i+1]*bet[i]) +       // T2(i+1,i)
                   Basics::abs(bet[i]*bet[i+1]) );                      // T2(i+2,i)
        normb = Basics::max(normb, nrm);
        i = k-2;
        nrm = sqrt(Basics::abs(bet[i-1]*bet[i-2]) +                     // T2(i-2,i)
                   Basics::abs(bet[i-1]*alp[i-1] + alp[i]*bet[i-1]) +   // T2(i-1,i)
                   bet[i-1]*bet[i-1] + alp[i]*alp[i] + bet[i]*bet[i] +  // T2(i,i)
                   Basics::abs(alp[i]*bet[i] + alp[i+1]*bet[i]) );      // T2(i+1,i)
        normb = Basics::max(normb, nrm);
        i = k-1;
        nrm = sqrt(Basics::abs(bet[i-1]*bet[i-2]) +                     // T2(i-2,i)
                   Basics::abs(bet[i-1]*alp[i-1] + alp[i]*bet[i-1]) +   // T2(i-1,i)
                   bet[i-1]*bet[i-1] + alp[i]*alp[i] );                 // T2(i,i)
        normb = Basics::max(normb, nrm);
    }
    return normb;
} */
