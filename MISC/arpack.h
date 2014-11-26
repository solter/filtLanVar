#ifndef ARPACK_H
#define ARPACK_H

#include <iostream>  // for cout, cerr, endl, etc.
#include <math.h>    // for sqrt, pow, etc.
#include "matkitdef.h"

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


struct ArpackOptions {
    Real tol;
    int numLanczosBasisVectors;

    // default constructor
    ArpackOptions () {
        tol = pow(Basics::machineEpsilon, 4.0/5.0);  // if it is not positive, it will be set to the default tolerance (machine epsilon) by dsaupd_() (for Real==double) or ssaupd_() (for Real==float)
        numLanczosBasisVectors = 0;  // the ncv in ARPACK; if 0, it will be set by 4*nev, with nev the number of desired eigenvalues
    }
};

struct ArpackInfo {
    int info;                              // error flag of ARPACK dsaupd_() routine
                                           // 0: normal exit
                                           // 1: maximum number of iterations taken
                                           // 2: deprecated, no longer an informational error
                                           // 3: no shifts could be applied during a cycle of the
                                           //    implicitly restarted Lanczos iteration;
                                           //    try increasing ncv
    int ierr;                              // error flag of ARPACK dseupd_() routine
                                           // 0: normal exit
    int numConvergedRitzValues;            // number of converged Ritz values to the requested tolerance
    int numImplicitRestarts;               // number of implicit restarts
    int numLanczosIter;                    // number of Lanczos iterations

    Real reverseCommunicationCpuTime;      // CPU time for reverse communication
    Real matrixVectorProductCpuTime;       // CPU time for matrix-vector products
    Real forEigenVectorsCpuTime;           // CPU time for obtaining eigenvectors
    int memoryForMatrixInBytes;            // memory required for storing the (sparse) matrix in bytes
    int memoryForArpackInBytes;            // memory required by ARPACK in bytes
};

// ArpackInfo EigenSolverARPACK(Vector &lambda, Matrix &U, const SymmetricMatrix &A, const int nev, const char which[]);
ArpackInfo ArpackEigenSolver(Vector &lambda, Matrix &U, const SparseMatrix &A, const int nev, const char which[], ArpackOptions &opts);
ArpackInfo ArpackEigenSolver(Vector &lambda, Matrix &U, const SparseMatrix &A, const int nev, const char which[]);  // default opts will be used

#endif
