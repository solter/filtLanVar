#ifndef SYMTRIEIG_H
#define SYMTRIEIG_H
/*! \file symtrieig.h
 *  Rountines for symmetric tridiagonal eigenvalue problems (for all eigenvalues and optionally all eigenvectors of a symmetric tridiagonal matrix).
 *
 *  The <b>LAPACK</b> routine <b>xSTEV</b> (i.e. <b>dstev</b>_() if <b>Real</b> is <b>double</b>, or <b>sstev</b>_() if <b>Real</b> is <b>float</b>) is required.
 *  The solvers declared here are indeed wrappers.
 */
// routines for symmetric eigenvalue problems, invoking LAPACK routine sstev_ (single precision) or dstev_ (double precision)

#include <iostream>  // for ostream
#include "matkitdef.h"

#ifdef USE_NAMESPACE
using namespace MATKIT;
namespace MATKIT {
#endif

class Vector;
class Matrix;

#ifdef USE_NAMESPACE
}  // end of namespace MATKIT
#endif



/** @defgroup group_symtrieig Eigenvalues / Eigenvectors of a Symmetric Tridiaognal Matrix
 *
 *  @{
 */

//! Compute all eigenvalues (but not eigenvectors) of a symmetric tridiagonal matrix <em>T</em>. This function invokes the <b>LAPACK</b> routine <b>dstev</b>_() (if <b>Real</b> is <b>double</b>) or <b>stev</b>_() (if <b>Real</b> is <b>float</b>).
/*! \param  n          is the dimension of the symmatrix tridiagonal matrix <em>T</em>.
 *  \param  diag[],sdiag[] define the symmetric tridiagonal matrix <em>T</em>, whose the diagonal elements are <em>diag</em>[<em>0,...,n-1</em>] and
 *          subdiagonal elements are <em>sdiag</em>[<em>0,...,n-2</em>]. In other words,
 *          <ul>
 *          <li> <em>T</em>(<em>i,i</em>)<em> = diag</em>[<em>i-1</em>]  for <em>i = 1,...,n</em>.
 *          <li> <em>T</em>(<em>i+1,i</em>)<em> = T</em>(<em>i,i+1</em>)<em> = sdiag</em>[<em>i-1</em>] for <em>i = 1,...n-1</em>.
 *          </ul>
 *  \param  eigVal     is the output vector of length <em>n</em> containing all eigenvalues of <em>T</em> in ascending order.
 *
 *  \return The flag returned by the <b>LAPACK</b> routine <b>dstev</b>_() (if <b>Real</b> is <b>double</b>) or <b>stev</b>_() (if <b>Real</b> is <b>float</b>).
 *
 *  \see    reportTroubleIfAny().
 */
// compute all eigenvalues (but not eigenvectors) of a symmetric tridiagonal matrix
// n is the dimension of the symmetric tridiagonal matrix
// diag[],sdiag[] define the symmetric tridiagonal matrix:
//     the diagonal elements are diag[0,...,n-1] in order and the subdiagonal elements are sdiag[0,...,n-2] in order
// eigVal is the output vector of length n containing all eigenvalues in ascending order
// the return value is the flag returned by the LAPACK routine dstev_() (if Real is double) or stev_() (if Real is float)
int SymmetricTridiagoanlEigenSolver(Vector &eigVal, mkIndex n, const Real *diag, const Real *sdiag);

//! Compute all eigenvalues and eigenvectors of a symmetric tridiagonal matrix <em>T</em>. This function invokes the <b>LAPACK</b> routine <b>dstev</b>_() (if <b>Real</b> is <b>double</b>) or <b>stev</b>_() (if <b>Real</b> is <b>float</b>).
/*! \param  n          is the dimension of the symmatrix tridiagonal matrix <em>T</em>.
 *  \param  diag[],sdiag[] define the symmetric tridiagonal matrix; the diagonal elements are <em>diag</em>[<em>0,...,n-1</em>] and
 *          subdiagonal elements are <em>sdiag</em>[<em>0,...,n-2</em>]. In other words,
 *          <ul>
 *          <li> <em>T</em>(<em>i,i</em>)<em> = diag</em>[<em>i-1</em>]  for <em>i = 1,...,n</em>.
 *          <li> <em>T</em>(<em>i+1,i</em>)<em> = T</em>(<em>i,i+1</em>)<em> = sdiag</em>[<em>i-1</em>] for <em>i = 1,...n-1</em>.
 *          </ul>
 *  \param  eigVal      is the output vector of length <em>n</em> containing all eigenvalues of <em>T</em> in ascending order.
 *  \param  eigVec      is the output <em>n</em>-by-<em>n</em> matrix with columns as eigenvectors of <em>T</em>, sorted with respect to <em>eigVal</em>.
 *
 *  \return The flag returned by <b>dstev</b>_() (if <b>Real</b> is <b>double</b>) or <b>stev</b>_() (if <b>Real</b> is <b>float</b>).
 *
 *  \see    reportTroubleIfAny().
 */
// compute all eigenvalues and eigenvectors of a symmetric tridiagonal matrix
// n is the dimension of the symmatrix tridiagonal matrix
// diag[],sdiag[] define the symmetric tridiagonal matrix:
//     the diagonal elements are diag[0,...,n-1] in order and the subdiagonal elements are sdiag[0,...,n-2] in order
// eigVal is the output vector of length n containing all eigenvalues in ascending order
// eigVec is the output n-by-n matrix with columns as eigenvectors, in the order as elements in eigVal
// the return value is the flag returned by the LAPACK routine dstev_() (if Real is double) or stev_() (if Real is float)
int SymmetricTridiagoanlEigenSolver(Vector &eigVal, Matrix &eigVec, mkIndex n, const Real *diag, const Real *sdiag);

//! Print out an error message to <b>std::ostream</b> <em>outerr</em> if <em>info\f$\neq\f$ 0</em> which signifies an error that occurred in the <b>LAPACK</b> routine <b>xSTEV</b>.
/*! \param  outerr      is an <b>std::ostream</b> for printing out the error message when <em>info\f$\neq\f$ 0</em>.
 *  \param  info        is the flag returned from the <b>LAPACK</b> routine <b>dstev</b>_() (if <b>Real</b> is <b>double</b>) or <b>sstev</b>_() (if <b>Real</b> is <b>float</b>).
 *          <ul>
 *          <li> =0:  successful exit  (nothing will be printed).
 *          <li> <0:  if <em>info=-i</em>, the <em>i</em>th argument had an illegal value.
 *          <li> >0:  if <em>info=i</em>, the algorithm failed to converge; <em>i</em> off-diagonal elements did not converge to zero.
 *          </ul>
 *  \param  n           is the dimension of the symmetric tridiagonal matrix.
 *          This information will also be printed out if <em>n</em> is provided (i.e. <em>n\f$\neq\f$ 0</em>).
 *
 *  \see    SymmetricTridiagoanlEigenSolver().
 */
// print out an error message to std::ostream outerr if info!=0 which signifies an error that occurred in the LAPACK routine xSTEV
// outerr is an std::ostream for printing out the error message when info!=0
// info is the flag returned from the LAPACK routine dstev_() (if Real is double) or sstev_() (if Real is float)
// = 0:  successful exit  (nothing will be printed)
// < 0:  if info = -i, the i-th argument had an illegal value
// > 0:  if info = i, the algorithm failed to converge; i off-diagonal elements did not converge to zero
// n is the dimension of the symmetric tridiagonal matrix; this information will also be printed out if n is provided (i.e. n!=0)
void reportTroubleIfAny(std::ostream &outerr, int info, mkIndex n=0);

/** @} */  // end of group_symtrieig


#endif
