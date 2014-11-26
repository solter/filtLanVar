#include <stdlib.h>  // for exit
#include <string.h>  // for memcpy
#include <time.h>    // for time_t, time, clock_t, clock, CLOCKS_PER_SEC, etc.
#include <iostream>  // for cout, cerr, endl, etc. (under namespace std)

#include "matkit.h"
#include "arpack.h"

using std::cout;
using std::endl;

#define FORTRAN_LOGICAL int
// with some compilers / machines, one may have to set Fortran logical as C long int, in which case the following line should be used instead of the above
//#define FORTRAN_LOGICAL long

#define USE_CSRMV
// A*v corresponds to CSCMV and v*A corresponds to CSRMV
// when A is (sparse) symmetric, then A*v equals v*A in exact arithmetic
// however CSRMV usually slightly faster than CSCMV
// "#define USE_CSRMV" invokes CSRMV

extern "C" {
    // routines from ARPACK
    void dsaupd_(int *ido,                 // 0 for first call to the reverse communication interface
                 char *bmat,               // 'I' for standard eigenvalue problem, or 'G' for generalized eigenvalue problem
                 int *n,                   // dimension of the eigenproblem
                 const char *which,        // size 2, "LA", "SA", "LM", "SM", or "BE"
                 const int *nev,           // number of eigenvalues wanted
                 double *tol,              // if less than or equal to 0, it will be set as DLAMCH('EPS') machine epsilon
                 double *resid,            // size *n, residual vector (probably from previous run), or a random vector (info = 0)
                 int *ncv,                 // # columns of V
                 double *V,                // V(*n,*ncv)
                 int *ldv,                 // leading dimension of V
                 int *iparam,              // size 11
                 int *ipntr,               // size 11
                 double *workd,            // size 3*(*n)
                 double *workl,            // size *lwork
                 int *lworkl,              // work size
                 int *info);

    void dseupd_(FORTRAN_LOGICAL *rvec,    // true for Ritz vectors, and false for Ritz values only
                 char *howmny,             // a character, 'A' for all *nev Ritz vectors, or
                                           // 'S' for some of the Ritz vectors, specified by the logial array select[]
                 FORTRAN_LOGICAL *select,  // size *ncv
                 double *D,                // size *nev, as output, D contains the Ritz values, in ascending order,
                                           // approximating the eigenvalues
                 double *Z,                // size *n, as output, Z contains the B-orthonormal Ritz vectors
                                           // if *rvec == false, Z is not referenced
                 int *ldz,                 // the leading dimension of Z
                                           // if Ritz vectors are desired, then *ldz >= max(1,*n)
                 double *sigma,            // a number represents the shift, if iparam[6]=3,4,5;
                                           // not referenced if iparam[6]==1,2 (iparam[6] is the model index)
                 // the following parameters are the same as those dsaupd_
                 char *bmat,               // 'I' for standard eigenvalue problem; 'G' for generalized eigenvalue problem
                 int *n,                   // dimension of the eigenproblem
                 const char *which,        // size 2, "LA", "SA", "LM", "SM", or "BE"
                 const int *nev,           // number of eigenvalues wanted
                 double *tol,              // if less than or equal to 0, it will be set as DLAMCH('EPS') machine epsilon
                 double *resid,            // size *n, residual vector (probably from previous run), or a random vector (info = 0)
                 int *ncv,                 // # columns of V
                 double *V,                // V(*n,*ncv)
                 int *ldv,                 // leading dimension of V
                 int *iparam,              // size 11
                 int *ipntr,               // size 11
                 double *workd,            // size 3*(*n)
                 double *workl,            // size *lwork
                 int *lworkl,              // work size
                 int *info);
}

ArpackInfo ArpackEigenSolver(Vector &lambda, Matrix &U, const SparseMatrix &A, const int nev, const char which[]) {
    ArpackOptions opts;
    return ArpackEigenSolver(lambda, U, A, nev, which, opts);
}

ArpackInfo ArpackEigenSolver(Vector &lambda, Matrix &U, const SparseMatrix &A, const int nev, const char which[], ArpackOptions &opts) {
    ArpackInfo outputInfo;
    int ido = 0;               // default has to be zero (initially)
    char bmat = 'I';           // 'I' for standard eigenvalue problem, or 'G' for generalized eigenvalue problem
    int n = A.Nrows();
    double *resid = new Real[n];
    int ncv = opts.numLanczosBasisVectors;
    if (ncv == 0)
        ncv = 2*nev;           // ncv is the largest number of basis vectors that will be used in the Implicitly Restarted Lanczos Process
    ncv = (ncv<n) ? ncv : n;
    int ldv = n;
    double *V = new Real[ldv*ncv];
    int iparam[11], ipntr[11];
    double *workd = new Real[3*n];
    int lworkl = ncv*(ncv+8);  // this is the value set in the ARPACK test driver dsdrv1.f
    double *workl = new Real[lworkl];
    outputInfo.info = 0;       // default is 0 for a random initial Lanczos vector
    outputInfo.memoryForMatrixInBytes = A.Nnz()*(sizeof(Real)+sizeof(int)) + (A.Ncols()+1)*sizeof(int);    // memory to store the matrix
    outputInfo.memoryForArpackInBytes = (n + ldv*ncv + 3*n + lworkl + 2*ncv)*sizeof(Real) + ncv*sizeof(FORTRAN_LOGICAL);  // memory for Lanczos iterations

    // specification of algorithm mode
    // the parameters follow the default values in the ARPACK test driver dsdrv1.f
    int ishfts = 1;
    int maxiter = 50*n;
    int model = 1;
    iparam[0] = ishfts;
    iparam[2] = maxiter;
    iparam[6] = model;

    clock_t reverseCommunicationClockTicks = 0;
    clock_t matrixVectorProductClockTicks = 0;
    for (int iter=1; iter<=maxiter; iter++) {
        clock_t start = clock();
        dsaupd_(&ido, &bmat, &n, which, &nev, &opts.tol, resid, &ncv, V, &ldv,
                iparam, ipntr, workd, workl, &lworkl, &outputInfo.info);
        reverseCommunicationClockTicks += (clock()-start);
        if (ido != -1 && ido != 1)
            break;
        start = clock();
        const Real *sv = workd+ipntr[0]-1;
        Real *sw = workd+ipntr[1]-1;
        #ifdef USE_CSRMV
            vector_spmatrix_multiplication((mkIndex)n, sv, A.Store(), A.RowIndex(), A.ColPointer(), sw);
        #else
            spmatrix_vector_multiplication((mkIndex)n, (mkIndex)n, A.Store(), A.RowIndex(), A.ColPointer(), sv, sw);
        #endif
        matrixVectorProductClockTicks += (clock()-start);
    }
    if (outputInfo.info < 0) {
        cout << "ArpackEigenSolver(Vector &, Matrix &, const SparseMatrix &, const int, const char []): fatal error, info=" << outputInfo.info << " from dsaupd_ !" << endl;
        cout << "Check the documentation in dsaupd_(). " << endl;
        exit(1);
    }
    // else no fatal error occurred

    // number of iterations, CPU time
    outputInfo.reverseCommunicationCpuTime = (Real)reverseCommunicationClockTicks / (Real)CLOCKS_PER_SEC;
    outputInfo.matrixVectorProductCpuTime = (Real)matrixVectorProductClockTicks / (Real)CLOCKS_PER_SEC;

    // compute eigenvalues and eigenvectors
    FORTRAN_LOGICAL rvec = 1;                            // want vectors (true)
    char howmny[] = "All";                               // for all nev Ritz vectors
    FORTRAN_LOGICAL *select = new FORTRAN_LOGICAL[ncv];  // allocate an array of size ncv, following ARPACK test driver dsdrv1.f
    Real sigma;                                          // shift for model 3,4,5; we use model==1 so not referenced
    Real *lam = new Real[2*ncv];                         // eigenvalues to be stored
    clock_t start = clock();
    dseupd_(&rvec, howmny, select, lam, V, &ldv, &sigma,
            &bmat, &n, which, &nev, &opts.tol, resid, &ncv, V, &ldv,
            iparam, ipntr, workd, workl, &lworkl, &outputInfo.ierr);
    outputInfo.forEigenVectorsCpuTime = (Real)(clock()-start) / (Real)CLOCKS_PER_SEC;
    outputInfo.numConvergedRitzValues = iparam[4];       // number of converged Ritz values to the requested tolerance
    opts.numLanczosBasisVectors = ncv;                   // number of Lanczos basis vectors
    outputInfo.numImplicitRestarts = iparam[2];          // number of implicit restarts
    outputInfo.numLanczosIter = iparam[8];               // number of iterations
    delete [] resid;
    delete [] workd;
    delete [] workl;
    delete [] select;
    int nevx = outputInfo.numConvergedRitzValues;
    if (nev < nevx || nevx < 0) {
       cout << "Warning: reported number of converged Ritz values is " << nevx << "!" << endl;
       nevx = nev;
    }

    int sz = n*nevx;
    Real *us = new Real[sz];
    Real *ss = us;
    Real *vs = V;
    int ii = nevx;
    while (ii--) {
        memcpy(ss, vs, n*sizeof(Real));
        ss += n;
        vs += ldv;
    }
    delete [] V;
    U.set(n, nevx, us);                                  // set eigenvectors
    Real *ls = new Real[nevx];
    memcpy(ls, lam, nevx*sizeof(Real));
    delete [] lam;
    lambda.set(nevx, ls);                                // set eigenvalues
    return outputInfo;
}
