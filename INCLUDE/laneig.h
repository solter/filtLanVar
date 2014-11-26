#ifndef LANEIG_H
#define LANEIG_H
/*! \file laneig.h
 *  Routines for the Lanczos method for extreme eigenvalues (and optionally the corresponding eigenvectors) of a symmetric matrix.
 */

#include <math.h>    // for sqrt, pow, etc.
#include <iostream>  // for ostream, cout, cerr, endl, etc. (under namespace std)

#ifdef USE_NAMESPACE
using namespace MATKIT;
namespace MATKIT {
#endif

class Vector;
class Matrix;
class SymmetricMatrix;
class SparseMatrix;

#ifdef USE_NAMESPACE
}  // end of namespace MATKIT
#endif



/** @defgroup group_laneig Lanczos Algorithm for Eigenvalues / Eigenvectors of a Symmetric Matrix
 *
 *  @{
 */

//! An instance of this class, taken by LanczosEigenSolver(), is a collection of options for the Lanczos procedure to solve symmetric eigenvalue problems.
// an instance of this class, taken by LanczosEigenSolver(), is a collection of options for the Lanczos procedure to solve symmetric eigenvalue problems
struct LanczosOptions {
    //! This parameter tells whether eigenvectors are wanted or not (default <b>true</b>).
    /*! Note that the standard Lanczos procedure computes the eigenvalues first.
     *  The computation of eigenvectors is unnecessary if only the eigenvalues are required.
     */
    // this parameter tells whether eigenvectors are wanted or not (default true)
    // note that the standard Lanczos procedure computes the eigenvalues first
    // the computation of eigenvectors is unnecessary if only the eigenvalues are required
    bool wantEigVec;
    //! Minimum number of Lanczos iterations.
    /*! This is the number of iterations performed before checking convergence.
     *
     *  The default value is <em>0</em>, which means that it will be determined by LanczosOptions::defaultMinIterFactor.
     */
    // minimum number of Lanczos iterations
    // this is the number of iterations performed before checking convergence
    // the default value is 0, which means that it will be determined by LanczosOptions::defaultMinIterFactor
    mkIndex minIter;
    //! Maximum number of Lanczos iterations.
    /*! This is the maximum number of iterations performed even if some desired eigenvalues do not yet converge.
     *
     *  The default value is <em>0</em>, which means that it will be determined by LanczosOptions::defaultMaxIterFactor.
     */
    // maximum number of Lanczos iterations
    // this is the maximum number of iterations performed even if some desired eigenvalues do not yet converge
    // the default value is 0, which means that it will be determined by LanczosOptions::defaultMaxIterFactor
    mkIndex maxIter;
    //! Extra number of Lanczos iterations (default <em>100</em>).
    /*! This is the number of Lanczos iterations performed after the Ritz values are `deemed' converged to the desired eigenvalues.
     *  There are two reasons to perform more Lanczos iterations.
     *  <ul>
     *  <li>  To improve the accuracy of the computed eigenvalues / eigenvectors.
     *  <li>  It is possible that there are missing eigenvalues. \n
     *  </ul>
     *  The second can be illustrated by the following example.
     *  When the <em>10</em> largest Ritz values are converged to 10 eigenvalues,
     *  it is possible that the <em>10</em> eigenvalues are not the largest,
     *  because one or more of the <em>10</em> largest eigenvalues do not yet show up as Ritz values.
     *  Nevertheless, these missing eigenvalues, if any, may come out by applying more iterations.
     *  This issue can be serious in a filtered Lanczos procedure with a high degree polynomial.
     *
     *  Finally, the phase of such extra iterations may be repeated for a few times until all desired eigenvalues are found.
     */
    // extraIter is the extra number of Lanczos iterations (default 100)
    // this is the number of Lanczos iterations performed after the Ritz values are `deemed' converged to the desired eigenvalues
    // there are two reasons to perform more Lanczos iterations
    // first, to improve the accuracy of the computed eigenvalues / eigenvectors
    // second, it is possible that there are missing eigenvalues
    // for instance, when the 10 largest Ritz values are converged to 10 eigenvalues,
    // it is possible that the 10 eigenvalues are not the largest,
    // because one or more of the 10 largest eigenvalues do not yet show up as Ritz values
    // nevertheless, these missing eigenvalues, if any, may come out by applying more iterations
    // this issue can be serious in a filtered Lanczos procedure with a high degree polynomial
    // finally, the phase of such extra iterations may be repeated for a few times until all desired eigenvalues are found
    mkIndex extraIter;
    //! The convergence is checked every `stride' iterations (default <em>10</em>).
    // the convergence is checked every `stride' iterations (default 10)
    mkIndex stride;
    //! The default LanczosOptions::minIter, if not set, is LanczosOptions::defaultMinIterFactor times the number of (largest or smallest) eigenvalues requested.
    /*! The default value of LanczosOptions::defaultMinIterFactor is <em>5</em>.
     */
    // The default LanczosOptions::minIter, if not set, is LanczosOptions::defaultMinIterFactor times the number of (largest or smallest) eigenvalues requested
    // the default value of LanczosOptions::defaultMinIterFactor is 5
    mkIndex defaultMinIterFactor;
    //! The default LanczosOptions::maxIter, if not set, is <em>500</em> plus LanczosOptions::defaultMaxIterFactor times the number of (largest or smallest) eigenvalues requested.
    /*! The default value of LanczosOptions::defaultMaxIterFactor is <em>50</em>.
     */
    // The default LanczosOptions::maxIter, if not set, is 500 plus LanczosOptions::defaultMaxIterFactor times the number of (largest or smallest) eigenvalues requested
    // the default value of LanczosOptions::defaultMaxIterFactor is 50
    mkIndex defaultMaxIterFactor;
    //! Tolerance of eigenvalues, defined in the relative and average sense, for convergence check (default \f$\sqrt[5]{\epsilon^4}\f$, with \f$\epsilon\f$ the machine epsilon).
    // tolerance of eigenvalues, defined in the relative and average sense, for convergence check (default eps^{0.8}, with eps the machine epsilon)
    Real tol;
    //! The initial Lanczos vector.
    /*! If <em>v0.Length</em>() is not the number of rows/columns of the matrix, then
     *  a random vector is generated for the initial Lanczos vector as in the default case.
     */
    // the initial Lanczos vector
    // if v0.Length() is not the number of rows/columns of the matrix, then
    // a random vector is generated for the initial Lanczos vector as the default case
    Vector v0;

    // variables for reorthogonalization
    //! A flag to tell which reorthogonalization scheme is used (default <em>1</em>). The value can be <em>0,1,2</em>.
    /*! <ul>
     *  <li>  <em>0</em> for no reorthogonalization (reserved option);
     *  <li>  <em>1</em> for partial reorthogonalization (default and recommended);
     *  <li>  <em>2</em> for full reorthogonalization.
     *  </ul>
     */
    // a flag to tell which reorthogonalization scheme is used (default 1):
    // 0 for no reorthogonalization (reserved option);
    // 1 for partial reorthogonalization (default and recommended);
    // 2 for full reorthogonalization
    int reorth;
    //! Reorthogonalization is doubled in the iterations with <em>nrm < doubleReorthGamma*nrm_old</em>, where <em>nrm_old</em> and <em>nrm</em> are the norms of the latest Lanczos vector before and after reorthogonalization.
    /*! The default value of <em>doubleReorthGamma</em> is \f$\sqrt{2}\f$.
     *  <ul>
     *  <li>  <em>doubleReorthGamma==0.0</em> (or smaller) means that reorthogonalization is never doubled.
     *  <li>  <em>doubleReorthGamma==1.0</em> (or larger) means that reorthogonalization, whenever it is performed, is always doubled.
     *  </ul>
     *
     *  \remark This parameter is effective only if LanczosOptions::reorth==<em>1</em> or LanczosOptions::reorth==<em>2</em> (i.e. partial or full reorthogonalization).
     */
    // reorthogonalization is doubled in the iterations with nrm < doubleReorthGamma*nrm_old
    // the default value of doubleReorthGamma is sqrt(0.5)
    // nrm_old and nrm are the norms of the latest Lanczos vector before and after reorthogonalization
    // doubleReorthGamma==0.0 (or smaller) means that reorthogonalization is never doubled
    // doubleReorthGamma==1.0 (or larger) means that reorthogonalization, whenever it is performed, is always doubled
    // this parameter is effective only if LanczosOptions::reorh==1 or LanczosOptions::reorth==2 (i.e. partial or full reorthogonalization)
    Real doubleReorthGamma;
    //! Local reorthogonalization is performed if <em>beta</em>[<em>j-1</em>]<<em>localReorthGamma*beta</em>[<em>j</em>] or <em>localReorthGamma>=1.0</em>, where <em>beta</em>[<em>j</em>] is the latest beta and <em>beta</em>[<em>j-1</em>] is the previous beta at <em>j</em>th Lanczos iteration.
    /*! The default value of <em>localReorthGamma</em> is sqrt(0.5). \n
     *  <em>localReorthGamma=0.0</em> (or smaller) means that local reorthogonalization is never performed.
     *
     *  \remark  This parameter is effective only if LanczosOptions::reorth==<em>1</em> (i.e. partial reorthogonalization) as input of LanczosEigenSolver().
     */
    // local reorthogonalization is performed if beta[j-1] < localReorthGamma*beta[j] or localReorthGamma >= 1.0,
    // the default value of localReorthGamma is sqrt(0.5)
    // where beta[j] is the latest beta and beta[j-1] is the previous beta at j-th Lanczos iteration
    // localReorthGamma=0.0 (or smaller) means that local reorthogonalization is never performed
    // this parameter is effective only if LanczosOptions::reorth==1 (i.e. partial reorthogonalization) as input of LanczosEigenSolver().
    Real localReorthGamma;
    //! A flag to tell which partial reorthogonalization strategy is used. The value can be <em>0,1,2,3</em>.
    /*! <ul>
     *  <li>  <em>0</em> for reorthogonalization against previous Lanczos vectors in groups; each group consists of
     *        <em>v</em>(<em>i</em>),<em>v</em>(<em>i+1</em>),...,<em>v</em>(<em>j</em>) with <em>omega</em>(<em>i</em>),<em>omega</em>(<em>i+1</em>),...,<em>omega</em>(<em>j</em>) all greater than <em>eta</em>, and
     *        there is <em>v</em>(<em>k</em>) with <em>omega</em>(<em>k</em>)><em>delta</em>, <em>i<=k<=j</em>.
     *  <li>  <em>1</em> for reorthogonalization against previous Lanczos vectors <em>v</em>(<em>i</em>) with <em>omega</em>(<em>i</em>)><em>eta</em>.
     *  <li>  <em>2</em> for reorthogonalization against previous Lanczos vectors in one group <em>v</em>(<em>i</em>),<em>v</em>(<em>i+1</em>),...,<em>v</em>(<em>j</em>) such that
     *        <em>i</em> is the smallest index for <em>omega</em>(<em>i</em>)><em>eta</em>, and <em>j</em> is the largest index for <em>omega</em>(<em>j</em>)><em>eta</em>.
     *  <li>  <em>3</em> for reorthogonalization against all previous Lanczos vectors.
     *  </ul>
     *
     *  <em>omega</em>(<em>i</em>) is the estimated orthogonal error of <em>v</em>(<em>i</em>) against the previous Lanczos vectors <em>v</em>(<em>1</em>),...,<em>v</em>(<em>i-1</em>).
     *  By default, the value of <em>eta</em> is <em>eps</em>^(<em>0.75</em>) and the value of <em>delta</em> is <em>eps</em>^(<em>0.5</em>),
     *  where <em>eps</em> is the machine epsilon.
     *
     *  \remark This parameter is effective only if LanczosOptions::reorth==<em>1</em> (i.e. partial reorthogonalization).
     */
    // a flag to tell which partial reorthogonalization strategy is used:
    // 0 for reorthogonalization against previous Lanczos vectors in groups; each group consists of
    //   v(i),v(i+1),...,v(j) with omega(i),omega(i+1),...,omega(j) all greater than eta, and
    //   there is v(k) with omega(k) > delta, i<=k<=j
    // 1 for reorthogonalization against previous Lanczos vectors v(i) with omega(i) > eta
    // 2 for reorthogonalization against previous Lanczos vectors in one group v(i),v(i+1),...,v(j) such that
    //   i is the smallest index for omega(i) > eta, and j is the largest index for omega(j) > eta
    // 3 for reorthogonalization against all previous Lanczos vectors
    // omega(i) is the estimated orthogonal error of v(i) against the previous Lanczos vectors v(1),...,v(i-1)
    // by default, the value of eta is eps^(0.75) and the value of delta is eps^(0.5),
    // where eps is the machine epsilon
    // this parameter is effective only if reorth==1 (i.e. partial reorthogonalization)
    int partialReorthStrategy;
    //! Set <em>checkReorth</em>=<b>true</b> for checking whether semi-orthogonality is preserved with partial reorthogonalization.
    /*! This is for verification purposes, and is expensive, comparable to full reorthogonalization.
     *  Therefore, <em>checkReorth</em>=<b>true</b> is not recommended for large problems.
     *
     *  By default, <em>checkReorth</em>==<b>false</b>.
     *
     *  \remark  This parameter is effective only with partial reorthogonalization (LanczosOptions::reorth==<em>1</em> as input of LanczosEigenSolver()).
     */
    // set checkReorth=true for checking whether semi-orthogonality is preserved with partial reorthogonalization
    // this is for verification purposes, and is expensive, comparable to full reorthogonalization
    // therefore, checkReorth=true is not recommended for large problems
    // by default, checkReorth==false
    // this parameter is effective only with partial reorthogonalization (LanczosOptions::reorth==1 as input of LanczosEigenSolver())
    bool checkReorth;
    //! This parameter determines how much the memory should be expanded when the allocated memory is not sufficient to store the Lanczos vectors and the elements alpha's and beta's which form the symmetric tridiagonal matrix.
    // memoryExpansionFactor determines how much the memory should be expanded
    // when the allocated memory is not sufficient to store the Lanczos vectors
    // and the elements alpha's and beta's which form the symmetric tridiagonal matrix
    Real memoryExpansionFactor;

    // cut points of eigenvalues
    //! Eigenvalues larger than <em>eigLowCut</em> in the lower end will not be computed.
    /*! This parameter is effective when <em>eigPart</em>[] is "SA" for smallest algebraic or when <em>eigPart</em>[] is "BE" for both ends.
     */
    // eigenvalues larger than eigLowCut in the lower end will not be computed
    // this parameter is effective when eigPart is "SA" for smallest algebraic or when eigPart is "BE" for both ends
    Real eigLowCut;
    //! Eigenvalues smaller than <em>eigHighCut</em> in the upper end will not be computed.
    /*! This parameter is effective when <em>eigPart</em>[] is "LA" for largest  algebraic or when <em>eigPart</em>[] is "BE" for both ends.
     */
    // eigenvalues smaller than eigHighCut in the upper end will not be computed
    // this parameter is effective when eigPart is "LA" for largest  algebraic or when eigPart is "BE" for both ends
    Real eigHighCut;
    //! A flag which indicates how to sort eigenvalues (default <em>-1</em>). The value should be no less than <em>-2</em>.
    /*! <ul>
     *  <li>  <em>-2</em> for `don't sort'.
     *  <li>  <em>-1</em> for the default setting:
     *        <ul>
     *        <li> If <em>eigPart</em>[] is "LA", <em>2</em> for to sort the eigenvalues in decreasing order.
     *        <li> If <em>eigPart</em>[] is "SA", <em>0</em> for to sort the eigenvalues in increasing order.
     *        <li> If <em>eigPart</em>[] is "BE", <em>0</em> for to sort the eigenvalues in increasing order.
     *        </ul>
     *  <li>  Otherwise eigSort should be >=<em>0</em>, in which case
     *        <ul>
     *        <li> The first  digit is <em>0</em> (or <em>1</em>) for to sort the eigenvalues in the lower end in increasing (or decreasing) order.
     *        <li> The second digit is <em>0</em> (or <em>1</em>) for to sort the eigenvalues in the higher end in increasing (or decreasing) order.
     *        <li> The third  digit is <em>0</em> (or <em>1</em>) for to list the eigenvalues in the lower end before (or after) those in the higher end.
     *        </ul>
     *  </ul>
     *  Here <em>eigPart</em>[] is the input string of LanczosEigenSolver().
     *
     *  \remark When <em>eigPart</em>[] is "LA", the first and the third digits do not matter.
     *  \remark When <em>eigPart</em>[] is "SA", the second and the third digits do not matter.
     *  \remark When <em>eigPart</em>[] is "BE", set <em>eigSort</em>=<em>0</em> (or <em>7</em>) for eigenvalues sorted in increasing (or decreasing) order.
     */
    // a flag which indicates how to sort eigenvalues (default -1):
    // -2 for `don't sort'
    // -1 for the default setting (2 for decreasing if eigPart[] is "LA"; 0 for increasing if eigPart[] is "SA"; 0 for increasing if eigPart[] is "BE")
    //  otherwise eigSort should be >= 0, in which case
    //        the first  digit is 0 (or 1) for to sort the eigenvalues in the lower  end in increasing (or decreasing) order
    //        the second digit is 0 (or 1) for to sort the eigenvalues in the higher end in increasing (or decreasing) order
    //        the third  digit is 0 (or 1) for to list the eigenvalues in the lower  end before (or after) those in the higher end
    //  Here eigPart[] is the input string of LanczosEigenSolver()
    //  note: when eigPart is "LA", the first and the third digits do not matter
    //        when eigPart is "SA", the second and the third digits do not matter
    //        when eigPart is "BE", set eigSort=0 (or 7) for eigenvalues sorted in increasing (or decreasing) order
    int eigSort;

    //! Diagnostic information display level (default <em>1</em>). The value can be <em>0,1,2</em>.
    /*! <ul>
     *  <li>  <em>0</em> means that no    information will be displayed.
     *  <li>  <em>1</em> means that basic information will be displayed.
     *  <li>  <em>2</em> means that more  information will be displayed.
     *  </ul>
     */
    // diagnostic information display level 0, 1 or 2 (default 1)
    // 0 means that no    information will be displayed
    // 1 means that basic information will be displayed
    // 2 means that more  information will be displayed
    int disp;
    //! A pointer to an output stream for diagnostic information. By default, it points to <em>std::cout</em>.
    // a pointer to an output stream for diagnostic information; by default, it points to std::cout
    std::ostream *out;
    //! A pointer to an output stream for error messages. By default, it points to <em>std::cerr</em>.
    // a pointer to an output stream for error messages; by default, it points to std::cerr
    std::ostream *err;

    //! A default constructor to set default options.
    // a default constructor to set default options
    LanczosOptions() {
        // default options
        wantEigVec = true;
        minIter = 0;      // 0 means that the default value will be set and used in LanczosEigenSolver
        maxIter = 0;      // 0 means that the default value will be set and used in LanczosEigenSolver
        extraIter = 100;  // extra number of iterations after convergence, to make sure no more eigenvalue to show up
        stride =  10;     // check convergence every `stride' iterations; if stride==0, it will be set as 10 in the routine LanczosEigenSolver

        defaultMinIterFactor = 5;
        defaultMaxIterFactor = 50;

        // reorthogonalization parameters
        reorth = 1;   // 1 for partial reorthogonalization; 2 for full reorthogonalization
        doubleReorthGamma = 1.0 / sqrt(2.0);
        localReorthGamma  = 1.0 / sqrt(2.0);
        partialReorthStrategy = 2;
        checkReorth = false;

        // memory expansion factor
        memoryExpansionFactor = 1.2;

        // tolerance for convergence check
        tol = pow(Basics::machineEpsilon, 4.0/5.0);

        // flag indicating how to sort the eigenvalues; see the description above
        eigSort = -1;

        // diagnostic information display parameters
        disp = 1;
        #ifdef USE_MEX
            out = &mexcout;
            err = &mexcout;
        #else
            out = &std::cout;
            err = &std::cerr;
        #endif

        // cut points of eigenvalues
        eigLowCut = Basics::infinity;    // the default means no cut
        eigHighCut = -Basics::infinity;  // the default means no cut
    }
};

//! LanczosEigenSolver() returns an instance of this class which gives the information of the Lanczos procedure.
// LanczosEigenSolver() returns an instance of this class which gives the information of the Lanczos procedure
struct LanczosInfo {
    //! Number of iterations performed in the Lanczos procedure.
    // number of iterations performed in the Lanczos procedure
    mkIndex numIter;
    //! CPU time (in total) for obtaining the next Krylov vector in each iteration, e.g. matrix-vector products.
    // CPU time (in total) for obtaining the next Krylov vector in each iteration, e.g. matrix-vector products
    Real forNextKrylovVectorCpuTime;
    //! CPU time (in total) for partial or full reorthogonalization.
    // CPU time (in total) for partial or full reorthogonalization
    Real reorthogonalizationCpuTime;
    //! CPU time (in total) for convergence check.
    // CPU time (in total) for convergence check
    Real convergenceCheckCpuTime;
    //! CPU time for obtaining eigenvectors, after the eigenvalues are obtained.
    // CPU time for obtaining eigenvectors, after the eigenvalues are obtained
    Real forEigenVectorsCpuTime;
    //! CPU time for the whole Lanczos procedure, including convergence check and getting eigenvectors (if requested).
    /*! \remark This is the total CPU time for computing the eigenvalues / eigenvectors,
     *          excluding reading the data matrix and reporting the result.
     */
    // CPU time for the whole Lanczos procedure,
    // including convergence check and getting eigenvectors (if requested);
    // this is the total CPU time for computing the eigenvalues / eigenvectors,
    // excluding reading the data matrix and reporting the result
    Real LanczosCpuTime;

    //! Maximum eigenvalue error.
    /*! The computation requires the eigenvectors, so this parameter is computed only if LanczosOptions::wantEigVec==<b>true</b> as the input of LanczosEigenSolver().
     */
    // maximum eigenvalue error
    // the computation requires the eigenvectors, so this parameter is computed only if LanczosOptions::wantEigVec==true as the input of LanczosEigenSolver()
    Real maxEigenvalueError;
    //! Maximum relative eigenvalue error.
    /*! The computation requires the eigenvectors, so this parameter is computed only if LanczosOptions::wantEigVec==<b>true</b> as the input of LanczosEigenSolver().
     */
    // maximum relative eigenvalue error
    // the computation requires the eigenvectors, so this parameter is computed only if LanczosOptions::wantEigVec==true as the input of LanczosEigenSolver()
    Real maxRelativeEigenvalueError;

    //! Memory (in bytes) required for the Lanczos iterations.
    // memory (in bytes) required for the Lanczos iterations
    mkIndex memoryForLanczosInBytes;

    //! <em>allEigenvaluesCheckedConverged</em> is <b>true</b> if all eigenvalues are checked converged.
    /*! Otherwise, <em>allEigenvaluesCheckedConverged</em> is <b>false</b>, which means that the maximum number of Lanczos iterations is reached.
     */
    // allEigenvaluesCheckedConverged is true if all eigenvalues are checked converged;
    // otherwise, allEigenvaluesCheckedConverged is false, which means that the maximum number of Lanczos iterations is reached
    bool allEigenvaluesCheckedConverged;

    // reorthogonalization cost information
    //! Number of iterations where reorthogonalization is performed.
    // number of iterations where reorthogonalization is performed
    mkIndex reorthIterCount;
    //! Total number of vectors against which the reorthogonalization is performed.
    // total number of vectors against which the reorthogonalization is performed
    mkIndex reorthVectorCount;
    //! The value is <em>reorthVectorCount</em> divided by the number of vectors reorthogonalized by full reorthogonalization without double reorthogonalization.
    /*! <ul>
     *  <li>  If no reorthogonalization (i.e. LanczosOptions::reorth==<em>0</em> as input of LanczosEigenSolver()) is performed, then <em>reorthVectorRate=0.0</em>.
     *  <li>  If full reorthogonalization (i.e. LanczosOptions::reorth==<em>2</em> as input of LanczosEigenSolver()) is performed without double reorthogonalization, then <em>reorthVectorRate=1.0</em>.
     *  <li>  If partial reorthogonalization (i.e. LanczosOptions::reorth==<em>1</em> as input of LanczosEigenSolver()) is performed, usually <em>0.0<reorthVectorRate<1.0</em>.
     *  </ul>
     */
    // the value is reorthVectorCount divided by the number of vectors reorthogonalized by full reorthogonalization without double reorthogonalization
    // if no reorthogonalization (i.e. LanczosOptions::reorth==0 as input of LanczosEigenSolver()) is performed, then reorthVectorRate = 0.0
    // if full reorthogonalization (i.e. LanczosOptions::reorth==2 as input of LanczosEigenSolver()) is performed without double reorthogonalization is performed, then reorthVectorRate = 1.0
    // if partial reorthogonalization (i.e. LanczosOptions::reorth==1 as input of LanczosEigenSolver()) is performed, usually 0.0 < reorthVectorRate < 1.0
    Real reorthVectorRate;
    //! Number of iterations where reorthogonalization is doubled.
    // number of iterations where reorthogonalization is doubled
    mkIndex localReorthIterCount;
    //! Number of iterations where local reorthogonalization is performed.
    // number of iterations where local reorthogonalization is performed
    mkIndex doubleReorthIterCount;

    //! Maximum ratio of orthogonality level.
    /*! Let <em>w0</em>(<em>i,j</em>) be the inner product of the <em>i</em>th and <em>j</em>th Lanczos vectors
     *  and <em>w0</em>(<em>j</em>) = max{<em>w0</em>(<em>i,j</em>)|<em>i<j</em>}. \n
     *  Also let w(j) be the estimated upper bound of w0(j). Then
     *  <em>maxOrthLevelRatio</em> = max{ <em>w0</em>(<em>j</em>)/<em>w</em>(<em>j</em>) | <em>1<=j<=numIter</em> }.
     *  \remark This parameter is computed only if the partial reorthogonalization is applied (i.e. LanczosOptions::reorth==<em>1</em> as input of LanczosEigenSolver()).
     */
    // maximum ratio of orthogonality level
    Real maxOrthLevelRatio;
    //! Minimum ratio of orthogonality level.
    /*! Let <em>w0</em>(<em>i,j</em>) be the inner product of the <em>i</em>th and <em>j</em>th Lanczos vectors
     *  and <em>w0</em>(<em>j</em>) = max{<em>w0</em>(<em>i,j</em>)|<em>i<j</em>}. \n
     *  Also let w(j) be the estimated upper bound of w0(j). Then
     *  <em>minOrthLevelRatio</em> = min{ <em>w0</em>(<em>j</em>)/<em>w</em>(<em>j</em>) | <em>1<=j<=numIter</em> }.
     *  \remark This parameter is computed only if the partial reorthogonalization is applied (i.e. LanczosOptions::reorth==<em>1</em> as input of LanczosEigenSolver()).
     */
    // minimum ratio of orthogonality level
    Real minOrthLevelRatio;
    // the above two parameters are computed only if the partial reorthogonalization is applied (i.e. LanczosOptions::reorth==1 as input of LanczosEigenSolver())
    // let w0(i,j) be the inner product of the i-th and j-th Lanczos vectors and w0(j) = max{w0(i,j)|i<j}
    // also let w(j) be the estimated upper bound of w0(j)
    // then minOrthLevelRatio = min{ w0(j)/w(j) | 1<=j<=numIter }, and
    //      maxOrthLevelRatio = max{ w0(j)/w(j) | 1<=j<=numIter }
};

//! The definition of a function / operator with input a vector and output a vector.
//! The operator should be self-adjoint if it is used for the Lanczos or conjugated gradient methods.
/*! <ul>
 *  <li>  In a standard Lanczos or conjugate gradient procedure, the input is <em>v</em> and output is <em>A*v</em>.
 *  <li>  In a filtered Lanczos or conjugate gradient procedure, the input is <em>v</em> and output is <em>p</em>(<em>A</em>)<em>*v</em>, with <em>p</em> a polynomial.
 *  </ul>
 *  In both cases, <em>A</em> is a symmetric matrix. It can be dense or sparse.
 */
// the definition of a function /operator with input a vector and output a vector
// the operator should be self-adjoint in the Lanczos methods
//     in a standard Lanczos or conjugate gradient procedure, the input is v and output is A*v;
//     in a filtered Lanczos or conjugate gradient procedure, the input is v and output is p(A)*v, with p(z) a polynomial,
// where A is a symmetric matrix; it can be dense or sparse
typedef Vector (*NEXT_VECTOR)(const Vector&);

//! Lanczos eigensolver, the most general form.
/*! \param  NextKrylovVector  is the function which defines a self-adjoint operator, e.g. <em>A*v</em> with <em>A</em> a symmetric matrix.
 *  \param  n         is the length of the input/output vectors of <em>NextKrylovVector</em>().
 *  \param  neigWanted  is the number of eigenvalues to be sought.
 *  \param  eigPart   is a string "LA", "SA", or "BE".
 *          <ul>
 *          <li> "LA" - largest algebraic, for largest eigenvalues;
 *          <li> "SA" - smallest algebraic, for smallest eigenvalues;
 *          <li> "BE" - both ends, one more from high end if nev is odd.
 *          </ul>
 *  \param  opts      is a collection of Lanczos options.
 *  \param  eigVal    is the output vector of length <em>neigWanted</em> containing the computed eigenvalues.
 *  \param  eigVec    is the output <em>n</em>-by-<em>neigWanted</em> matrix with columns as eigenvectors, in the order as elements in <em>eigVal</em>, if <em>opts.wantEigVec</em>==<b>true</b>.
 *
 *  \return An instance of LanczosInfo which gives the information of the Lanczos procedure.
 */
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
LanczosInfo LanczosEigenSolver(Vector &eigVal, Matrix &eigVec, NEXT_VECTOR NextKrylovVector, mkIndex n, mkIndex neigWanted, const char eigPart[], LanczosOptions &opts);
//! The same as LanczosEigenSolver(), but with the default LanczosOptions.
// the same as LanczosEigenSolver(), but with the default LanczosOptions
LanczosInfo LanczosEigenSolver(Vector &eigVal, Matrix &eigVec, NEXT_VECTOR NextKrylovVector, mkIndex n, mkIndex neigWanted, const char eigPart[]);
//! The same as LanczosEigenSolver(), but with the <em>NextKrylovVector</em>() defined as <em>A*v</em>.
// the same as LanczosEigenSolver(), but with the NextKrylovVector() defined as A*v
LanczosInfo LanczosEigenSolver(Vector &eigVal, Matrix &eigVec, const SymmetricMatrix &A, mkIndex neigWanted, const char eigPart[], LanczosOptions &opts);
//! The same as LanczosEigenSolver(), but with the <em>NextKrylovVector</em>() defined as <em>A*v</em> and with default LanczosOptions.
// the same as LanczosEigenSolver(), but with the NextKrylovVector() defined as A*v and with the default LanczosOptions
LanczosInfo LanczosEigenSolver(Vector &eigVal, Matrix &eigVec, const SymmetricMatrix &A, mkIndex neigWanted, const char eigPart[]);
//! The same as LanczosEigenSolver(), but with the <em>NextKrylovVector</em>() defined as <em>A*v</em>.
// the same as LanczosEigenSolver(), but with the NextKrylovVector() defined as A*v
LanczosInfo LanczosEigenSolver(Vector &eigVal, Matrix &eigVec, const SparseMatrix &A, mkIndex neigWanted, const char eigPart[], LanczosOptions &opts);
//! The same as LanczosEigenSolver(), but with the <em>NextKrylovVector</em>() defined as <em>A*v</em> and with default LanczosOptions.
// the same as LanczosEigenSolver(), but with the NextKrylovVector() defined as A*v and with the default LanczosOptions
LanczosInfo LanczosEigenSolver(Vector &eigVal, Matrix &eigVec, const SparseMatrix &A, mkIndex neigWanted, const char eigPart[]);

/** @} */  // end of group_laneig


#endif  // end of #ifdef LANEIG_H
