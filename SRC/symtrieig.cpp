//============================================================================
// Routines for computing eigenvalues of a symmetric tridiagonal matrix.
// They are wrappers of the LAPACK routine dstev_() or sstev_()
//============================================================================

#include <string.h>  // for memcpy, strcmp, strncmp, etc.
#include <iostream>  // for cout, cerr, endl, etc. (under namespace std)

#include "matkit.h"
#include "symtrieig.h"

#ifdef USE_MWLAPACK
    #include "lapack.h"
#else
    extern "C" {
        // interface of the LAPACK routines "sstev" (single precision) and "dstev" (double precision) for
        // all eigenvalues (and optional eigenvectors) of a symmetric tridiagonal matrix
        #ifdef USE_SINGLE
            void sstev_(char *jobz, int *n, float  *diagonal, float  *subdiagonal, float  *V, int *ldz, float  *work, int *info);
        #else
            void dstev_(char *jobz, int *n, double *diagonal, double *subdiagonal, double *V, int *ldz, double *work, int *info);
        #endif
    }
#endif

using std::endl;


// compute all eigenvalues (but not eigenvectors) of a symmetric tridiagonal matrix
// n is the dimension of the symmetric tridiagonal matrix
// diag[],sdiag[] define the symmetric tridiagonal matrix:
//     the diagonal elements are diag[0,...,n-1] in order and the subdiagonal elements are sdiag[0,...,n-2] in order
// eigVal is the output vector of length n containing all eigenvalues in ascending order
// the return value is the flag returned by the LAPACK routine dstev_() (if Real is double) or stev_() (if Real is float)
int SymmetricTridiagoanlEigenSolver(Vector &eigVal, mkIndex n, const Real *diag, const Real *sdiag) {
    char jobz = 'N';  // compute eigenvalues only
    #ifdef USE_MWLAPACK
        mwSignedIndex nn = (mwSignedIndex)n;
        mwSignedIndex ldz = nn;
        mwSignedIndex info;
    #else
        int nn = (int)n;
        int ldz = n;
        int info;  // output flag
    #endif

    // copy diagonal and subdiagonal elements to alp and bet
    Real *alp = new Real[n];
    Real *bet = new Real[n-1];
    memcpy(alp, diag, n*sizeof(Real));
    memcpy(bet, sdiag, (n-1)*sizeof(Real));

    // allocate storage for computation
    Real *sv = NULL;    // requires no additional storage since eigenvectors are not wanted
    Real *work = NULL;  // requires no additional storage since eigenvectors are not wanted
    #ifdef USE_SINGLE
        // single precision
        #ifdef USE_MWLAPACK
            sstev (&jobz, &nn, alp, bet, sv, &ldz, work, &info);
        #else
            sstev_(&jobz, &nn, alp, bet, sv, &ldz, work, &info);
        #endif
    #else
        // double precision
        #ifdef USE_MWLAPACK
            dstev (&jobz, &nn, alp, bet, sv, &ldz, work, &info);
        #else
            dstev_(&jobz, &nn, alp, bet, sv, &ldz, work, &info);
        #endif
    #endif

    // free memory
    delete [] bet;
    delete [] work;

    // set output
    eigVal.set(n, alp);

    // return info
    #ifdef USE_MWLAPACK
        return (int)info;
    #else
        return info;
    #endif
}

// compute all eigenvalues and eigenvectors of a symmetric tridiagonal matrix
// n is the dimension of the symmetric tridiagonal matrix
// diag[],sdiag[] define the symmetric tridiagonal matrix:
//     the diagonal elements are diag[0,...,n-1] in order and the subdiagonal elements are sdiag[0,...,n-2] in order
// eigVal is the output vector of length n containing all eigenvalues in ascending order
// eigVec is the output n-by-n matrix with columns as eigenvectors, in the order as elements in eigVal
// the return value is the flag returned by the LAPACK routine dstev_() (if Real is double) or stev_() (if Real is float)
int SymmetricTridiagoanlEigenSolver(Vector &eigVal, Matrix &eigVec, mkIndex n, const Real *diag, const Real *sdiag) {
    char jobz = 'V';  // compute eigenvalues and eigenvectors
    #ifdef USE_MWLAPACK
        mwSignedIndex nn = (mwSignedIndex)n;
        mwSignedIndex ldz = nn;
        mwSignedIndex info;
    #else
        int nn = (mkIndex)n;
        int ldz = n;
        int info;  // output flag
    #endif

    // copy diagonal and subdiagonal elements to alp and bet
    Real *alp = new Real[n];
    Real *bet = new Real[n-1];
    memcpy(alp, diag, n*sizeof(Real));
    memcpy(bet, sdiag, (n-1)*sizeof(Real));

    // allocate storage for computation
    Real *sv = new Real[n*n];
    Real *work = new Real[2*n-2];

    #ifdef USE_SINGLE
        // single precision
        #ifdef USE_MWLAPACK
            sstev (&jobz, &nn, alp, bet, sv, &ldz, work, &info);
        #else
            sstev_(&jobz, &nn, alp, bet, sv, &ldz, work, &info);
        #endif
    #else
        // double precision
        #ifdef USE_MWLAPACK
            dstev (&jobz, &nn, alp, bet, sv, &ldz, work, &info);
        #else
            dstev_(&jobz, &nn, alp, bet, sv, &ldz, work, &info);
        #endif
    #endif

    // free memory
    delete [] bet;
    delete [] work;

    // set output
    eigVal.set(n, alp);
    eigVec.set(n, n, sv);

    // return info
    #ifdef USE_MWLAPACK
        return (int)info;
    #else
        return info;
    #endif
}

// print out an error message to ostream outerr if info!=0 which signifies an error that occurred in the LAPACK routine xSTEV
// outerr is an ostream for printing out the error message when info!=0
// info is the flag returned from the LAPACK routine dstev_() (if Real is double) or sstev_() (if Real is float)
// = 0:  successful exit  (nothing will be printed)
// < 0:  if info = -i, the i-th argument had an illegal value
// > 0:  if info = i, the algorithm failed to converge; i off-diagonal elements did not converge to zero
// n is the dimension of the symmetric tridiagonal matrix; this information will also be printed out if n is provided (i.e. n!=0)
void reportTroubleIfAny(std::ostream &outerr, int info, mkIndex n) {
    if (info < 0) {
        // not supposed to happen; arguments should have been properly passed
        outerr << "Error: the " << -info;
        if (info == -1)
            outerr << "st" << endl;
        else if (info == -2)
            outerr << "nd" << endl;
        else if (info == -3)
            outerr << "rd" << endl;
        else
            outerr << "th" << endl;
        #ifdef USE_SINGLE
            outerr << " argument of sstev_ had an illegal value!" << endl;
        #else
            outerr << " argument of dstev_ had an illegal value!" << endl;
        #endif
    }
    else if (info > 0) {
        #ifdef USE_SINGLE
            outerr << "Warning: sstev_ failed to converge; " << info << " off-diagonal elements";
        #else
            outerr << "Warning: dstev_ failed to converge; " << info << " off-diagonal elements";
        #endif
        if (n != 0u)
            outerr << " (matrix dimension: " << n << ")";
        outerr << " did not converge to zero!" << endl;
    }
}
