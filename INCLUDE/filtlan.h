/*! \mainpage FILTLAN
 *
 *
 *  \section secwhat What is FILTLAN?
 *
 *  The <b>FILTLAN</b> library is a collection of <b>C</b>/<b>C++</b> routines of polynomial filtered Lanczos methods
 *  for computing extreme and interior eigenvalues of a symmetric matrix.
 *  <b>FILTLAN</b> does not handle generalized or other non-standard eigenvalue problems.
 *
 *  The ingredients include:
 *  <ul>
 *  <li> The Lanczos algorithm for symmetric eigenvalue computations.
 *  <li> A polynomial filter to transform the spectrum for the fast convergence of the desired eigenvalues. \n
 *       This polynomial is obtained from applying a conjugate-residual-type algorithm in polynomial space.
 *  <li> Partial reorthogonalization for maintaining semi-orthogonality among the Lanczos vectors.
 *  </ul>
 *
 *  Note:
 *  <ul>
 *  <li> This library uses <b>MATKIT</b> for basic matrix computations.
 *  <li> It also needs <b>LAPACK</b> for the routines <b>dstev</b>_() and <b>sstev</b>_() for eigenvalues of a symmetric tridiagonal matrix.
 *  </ul>
 *
 *  \section secfeedback Your feedback is very welcome!
 *
 *  Questions, comments, and complaints, please send us email to {<tt>hrfang</tt>, <tt>saad</tt>} (at) <tt>cs</tt> dot <tt>umn</tt> dot <tt>edu</tt>.
 *
 *  A technical report desribing the techniques used in the package can be found here: 
 *
 *  ``<em>A Filtered Lanczos Procedure for Extreme and Interior Eigenvalue Problems,</em>'' \n
 *  Haw-ren Fang and Yousef Saad,
 *  University of Minnesota Technical Report, 2011.
 */
#ifndef FILTLAN_H
#define FILTLAN_H
/*! \file filtlan.h
 *  Routines for the Lanczos method with polynomial filtering
 *  for extreme / interior eigenvalues and the corresponding eigenvectors of a symmetric matrix.
 */

#include "laneig.h"
#include "polyfilt.h"

#ifdef USE_NAMESPACE
using namespace MATKIT;
namespace MATKIT {
#endif

class Vector;
class Matrix;
class SparseMatrix;

#ifdef USE_NAMESPACE
}  // end of namespace MATKIT
#endif


/** @defgroup group_filtlan Filtered Lanczos Algorithm for Eigenvalues / Eigenvectors of a Symmetric Matrix
 *
 *  @{
 */

//! An instance of this class, taken by FilteredLanczosEigenSolver(), is a collection of options for the filtered Lanczos procedure to solve symmetric eigenvalue problems.
// an instance of this class, taken by FilteredLanczosEigenSolver(), is a collection of options for the filtered Lanczos procedure to solve symmetric eigenvalue problems
struct FilteredLanczosOptions: public LanczosOptions {
    //! This parameter tells whether eigenvalues are wanted or not (default <b>true</b>).
    /*! Note that the filtered Lanczos procedure computes the eigenvectors first.
     *  The computation of eigenvalues is unnecessary if only the invariant subspace is needed.
     *  \remark  The option <em>wantEigVec</em> inherited from LanczosOptions is no longer effective in the filtered Lanczos procedure.
     */
    // this parameter tells whether eigenvalues are wanted or not (default true)
    // note that the filtered Lanczos procedure computes the eigenvectors first
    // the computation of eigenvalues is unnecessary if only the invariant subspace is required
    // note that the option wantEigVec inherited from LanczosOptions is no longer effective in the filtered Lanczos procedure
    bool wantEigVal;
    //! Number of eigenvalues/eigenvectors desired.
    /*! If <em>neigWanted==0</em> is given (the default), then it means that all eigenvalues/eigenvectors in the requested interval are to be sought.
     *  Otherwise, let <em>nev</em> be the number of eigenvalues in the requested interval. Then
     *  <ul>
     *  <li>  If <em>nev <= neigWanted</em>, then the iteration continues until all <em>nev</em> eigenvalues/eigenvectors in the requested interval converge.
     *  <li>  If <em>nev > neigWanted</em>, then the iteration stops when <em>neigWanted</em> eigenvalues in the requested interval converge.
     *  </ul>
     */
    // neigWanted is the number of eigenvalues/eigenvectors desired
    // if nev == 0 is given (the default), then it means that all eigenvalues/eigenvectors in the requested interval are to be sought
    // otherwise, let nev be the number of eigenvalues in the requested interval
    //    if nev <= neigWanted, then the iteration continues until all nev eigenvalues/eigenvectors in the requested interval converge
    //    if nev >  neigWanted, then the iteration stops when neigWanted eigenvalues in the requested interval converge
    mkIndex neigWanted;
    //! Number of Lanczos iterations to determine an interval which (tightly) contains all eigenvalues (default <em>30</em>).
    /*! The filtered Lanczos procedure requires an interval which (tightly) contains all eigenvalues.
     *  If it is not provided, then we use a standard Lanczos procedure with <em>numIterForEigenRange</em> iterations to compute one.
     *  The higher the <em>numIterForEigenRange</em>, the tighter the interval.
     */
    // number of Lanczos iterations to determine an interval which (tightly) contains all eigenvalues (default 30)
    // the filtered Lanczos procedure requires an interval which (tightly) contains all eigenvalues
    // if it is not provided, then we use a standard Lanczos procedure with numIterForEigenRange iterations to compute one
    // the higher the numIterForEigenRange, the tighter the interval
    mkIndex numIterForEigenRange;

    //! A collection of options to determine the intervals.
    // a collection of options to determine the intervals
    IntervalOptions intervalOpts;

    //! A constructor to set default options.
    // a constructor to set default options
    FilteredLanczosOptions() {
        wantEigVal = true;          // filtered Lanczos computes eigenvectors first;
                                    // if only the invariant subspace is wanted, then the computation of eigenvalues is not required,
                                    // in which case wantEigVal can be set false
        neigWanted = 0;             // means to compute all eigenvalues in the given range
                                    // if neigWanted > 0, then quit the Lanczos loop when neigWanted eigenvalues are found in the specified range,
                                    //                    or when all ( < neigWanted ) eigenvalues are found in the specified range
        numIterForEigenRange = 30;  // number of Lanczos iterations to determine the eigen-range
        extraIter = 20;             // laneig has the default extraIter==60
                                    // filtlan requires less extraIter for comparable accuracy of eigenvalues
        defaultMinIterFactor = 2;   // likewise, filtlan can have less defaultMinIterFactor,
        defaultMaxIterFactor = 30;  // and also defaultMaxIterFactor than laneig
    }
};

//! FilteredLanczosEigenSolver() returns an instance of this class which gives the information of the filtered Lanczos procedure.
// FilteredLanczosEigenSolver() returns an instance of this class which gives the information of the filtered Lanczos procedure
struct FilteredLanczosInfo: public LanczosInfo {
    //! Information about the Lanczos procedure to determine the range of eigenvalues by approximating the largest eigenvalue and the smallest eigenvalue, if this range is not provided.
    // information about the Lanczos procedure to determine the range of eigenvalues by approximating the largest eigenvalue and the smallest eigenvalue, if this range is not provided
    LanczosInfo forEigenRangeInfo;
    //! CPU time (in seconds) to determine the range of eigenvalues, if this range is not provided.
    // CPU time (in seconds) to determine the range of eigenvalues, if this range is not provided
    Real forEigenRangeCpuTime;
    //! The range of eigenvalues is partitioned into intervals which determine the base filter approximated by a polynomial filter.
    /*! The <em>j</em>th interval is [<em>intervals</em>(<em>j</em>),<em>intervals</em>(<em>j+1</em>)).
     */
    // the range of eigenvalues is partitioned into intervals which determine the base filter approximated by a polynomial filter
    // the j-th interval is [intervals(j),intervals(j+1))
    Vector intervals;
};


//! This routine computes the eigenvalues in the requested window and the corresponding eigenvectors.
/*! <ul>
 *  <li>  The window can contain the largest or smallest eigenvalues, in which case the requested are the extreme eigenvalues.
 *  <li>  The window can also be well inside the interval of the spectrum, in which case it is often referred to as an interior eigenvalue problem.
 *  </ul>
 *
 *  \param  eigVal   is a vector which contains the computed eigenvalues (if <em>opts.wantEigVal</em>==<b>true</b>).
 *  \param  eigVec   is a matrix whose columns are computed eigenvectors. \n
 *          More precisely, <em>eigVec</em>(:,<em>j</em>) is the eigenvector corresponding to the eigenvalue <em>eigVal</em>(<em>j</em>).
 *  \param  S        is the sparse matrix whose eigenvalues / eigenvectors are to be sought.
 *  \param  frame    is a vector of 2 or 4 elements:
 *          <ul>
 *          <li>  If <em>frame</em> consists of 2 elements, then [<em>frame</em>(<em>1</em>),<em>frame</em>(<em>2</em>)] is the window in which the eigenvalues are sought.
 *          <li>  If <em>frame</em> consists of 4 elements, then [<em>frame</em>(<em>2</em>),<em>frame</em>(<em>3</em>)] is the window in which the eigenvalues are sought, and
 *                   [<em>frame</em>(<em>1</em>),<em>frame</em>(<em>4</em>)] is the interval which (tightly) contains all eigenvalues.
 *          </ul>
 *          In the former case, an interval which contains (tightly) all eigenvalues is not given, and will be determined by a standard Lanczos procedure.
 *          See also FilteredLanczosOptions::numIterForEigenRange.
 *  \param  polyDegree  is the degree of <em>s</em>(<em>z</em>), with <em>z*s</em>(<em>z</em>) the polynomial filter.
 *  \param  baseDegree  is the left-and-right degree of the base filter for each interval.
 *  \param  opts     is a collection of options for the filtered Lanczos procedure.
 *
 *  \return An instance of class FilteredLanczosInfo which gives the information of the filtered Lanczos procedure.
 */
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
                                       mkIndex polyDegree, mkIndex baseDegree, FilteredLanczosOptions &opts);

/** @} */  // end of group_filtlan


#endif
