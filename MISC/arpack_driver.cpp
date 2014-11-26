//============================================================================
// This arpack driver for computing extreme eigenvalues of a symmetric matrix
// using the arpack routines DSAUPD and DSEUPD.
//
// coded by H.-r. Fang, last update May, 2012.
//============================================================================

#include <stdlib.h>  // for exit, atof, atoi, atol, srand, rand, etc.
#include <string.h>  // for memcpy, strcmp, strncmp, etc.
#include <time.h>    // for time_t, time, clock_t, clock, CLOCKS_PER_SEC, etc.
#include <iostream>  // for cout, cerr, endl, etc. (under namespace std)
#include <iomanip>   // for setw, setprecision, etc. (under namespace std)

#include "matkit.h"
#include "arpack.h"

using std::cout;
using std::cerr;
using std::endl;

using std::setw;
using std::setprecision;
using std::scientific;

#ifdef USE_NAMESPACE
using namespace MATKIT;
#endif

void printUsageAndExit(char* cmd);


int main(int argc, char *argv[]) {
    char *cmd = argv[0];
    if (argc <= 2                                  // argc<=2 means no input matrix is specified
    || !strcmp(argv[1], "-?") || !strcmp(argv[1], "-h") || !strcmp(argv[1], "-help") || !strcmp(argv[1], "--help"))
        printUsageAndExit(cmd);                    // print usage and exit

    // declare variables for input parameters
    bool printEigVal = true;                       // whether to print eigenvalues or not
    bool printEigVec = false;                      // whether to print eigenvectors or not
    ArpackOptions opts;                            // a few ARPACK options
    mkIndex neigWanted = 6;                        // number of eigenvalues
    char defaultEigPart[] = "LM";                  // default part of desired eigenvalues
    char *eigPart = defaultEigPart;                // a string specifying which part of eigenvalues is desired

    // parse input parameters
    for (int i=2; i<argc; i++) {
        if (!strcmp(argv[i], "-srand"))
            srand(time(NULL));                              // set random seed as time(NULL)
        else if (!strncmp(argv[i], "-srand=", 7))
            srand(atoi(argv[i]+7));                         // random seed (default 1)
        else if (!strncmp(argv[i], "-printeig=", 10)) {
            if (atoi(argv[i]+10))                           // whether to print eigenvalue or not
                printEigVal = true;
            else
                printEigVal = false;
            if (atoi(argv[i]+10) >= 2)                      // whether to print eigenvector or not
                printEigVec = true;
            else
                printEigVec = false;
        }
        else if (!strncmp(argv[i], "-nev=", 5))
            neigWanted = atoi(argv[i]+5);                   // number of eigenvalues desired (default 6)
        else if (!strncmp(argv[i], "-part=", 6))
            eigPart = argv[i]+6;                            // this string specifies which part of eigenvalues are desired
                                                            // "LA", largest algebraic, for largest eigenvalues (default)
                                                            // "SA", smallest algebraic, for smallest eigenvalues
                                                            // "BE", both ends, one more from high end if nev is odd
        else if (!strncmp(argv[i], "-ncv=", 5))
            opts.numLanczosBasisVectors = atoi(argv[i]+5);  // maximum number of Lanczos basis vectors (default 2*nev)
        else if (!strncmp(argv[i], "-tol=", 5))
            opts.tol = atof(argv[i]+5);                     // tolerance for convergence test
        else {
            cerr << "Error: the parameter \"" << argv[i] << "\" is not recognized!" << endl;
            cerr << "For more information, use \"" << cmd << " --help\"." << endl;
            exit(1);
        }
    }

    if (neigWanted == 0) {
        cerr << "ARPACK driver: the number of desired eigenvalues is not specified!" << endl;
        exit(1);
    }

    // open the matrix file and read the matrix
    char *fileName = argv[1];
    cout << "input matrix file: " << fileName << endl;
    SparseMatrix A = SparseMatrix::mmread(fileName);
    cout << "n=" << A.Nrows() << ", nnz=" << A.Nnz() << endl;

    // print some input information
    cout << "ARPACK symmetric eigenvalue solver:" << endl;
    cout << "nev=" << neigWanted << ", eigPart=" << eigPart << endl;

    // compute eigenvalues by ARPACK
    clock_t start = clock();
    time_t start2 = time(NULL);
    Vector lambda;
    Matrix V;
    ArpackInfo outputInfo = ArpackEigenSolver(lambda, V, A, neigWanted, eigPart, opts);
    Real totalCpuTime = (Real)(clock()-start) / (Real)CLOCKS_PER_SEC;
    time_t wallClockTime = time(NULL) - start2;

    // report the result
    cout << "========================================================" << endl;
    if (outputInfo.info == 1) {
        cout << "maximum number of iterations taken (info==1)" << endl;
    }
    else if (outputInfo.info == 3) {
        cout << "no shifts could be applied during a cycle of the implicitly restarted Lanczos iteration (info==3)" << endl;
        cout << "try increasing ncv" << endl;
    }
    if (outputInfo.ierr < 0) {
        cout << "Error on output of ARPACK dseupd_ routine (" << outputInfo.ierr << ")" << endl;
    }
    cout << "number of converged Ritz values: " << outputInfo.numConvergedRitzValues << endl;
    cout << "number of Lanczos basis vectors: " << opts.numLanczosBasisVectors << endl;
    cout << "number of Lanczos iterations: " << outputInfo.numLanczosIter << endl;
    cout << "number of implicit restarts: " << outputInfo.numImplicitRestarts << endl;
    cout << "eigenvalue tolerance: " << opts.tol << endl;
    cout << "memory for storing matrix: " << outputInfo.memoryForMatrixInBytes << " bytes (" << Real(outputInfo.memoryForMatrixInBytes)/(Real)(1024*1024) << " MB)" << endl;
    cout << "memory required by ARPACK: " << outputInfo.memoryForArpackInBytes << " bytes (" << Real(outputInfo.memoryForArpackInBytes)/(Real)(1024*1024) << " MB)" << endl;
    cout << "CPU time for reverse communication: " << outputInfo.reverseCommunicationCpuTime << " seconds" << endl;
    cout << "CPU time for matrix-vector products: " << outputInfo.matrixVectorProductCpuTime << " seconds" << endl;
    cout << "CPU time used for obtaining eigenvectors: " << outputInfo.forEigenVectorsCpuTime << " seconds" << endl;
    cout << "total CPU time used: " << totalCpuTime << " seconds" << endl;
    cout << "wall clock time: " << wallClockTime << " seconds" << endl;
    cout << lambda.Length() << " eigenvalues found." << endl;
    Matrix E = A*V-V*lambda.spdiag();
    cout << "||A*V-V*S||_F=" << E.normFrobenius() << endl;
    cout << "||V'*V-I||_F=" << (V.transpose()*V-SparseMatrix::eye(lambda.Length())).normFrobenius() << endl;
    #ifdef USE_SINGLE
        const unsigned nprec = 7, nwidth= 17;
    #else
        const unsigned nprec = 15, nwidth= 25;
    #endif
    if (printEigVal) {
        cout << "========================================================" << endl;
        cout << scientific << setprecision(nprec);
        cout << "eig id          eigenvalue             residual norm" << endl;
        for (mkIndex i=1; i<=lambda.Length(); i++) {
            cout << setw(6) << i << setw(nwidth) << lambda(i) << setw(nwidth) << E.column(i).norm2() << endl;
        }
    }
    if (printEigVec) {
        cout << "========================================================" << endl;
        cout << scientific << setprecision(nprec);
        cout << "eig id          eigenvector" << endl;
        for (mkIndex i=1; i<=lambda.Length(); i++) {
            cout << setw(6) << i;
            for (mkIndex j=1; j<=V.Nrows(); j++)
                cout << setw(nwidth) << V(j,i);
            cout << endl;
        }
    }

    return 0;
}

void printUsageAndExit(char* cmd) {
    cout << "Usage: " << cmd << " MATRIX_FILE [OPTION] [OPTION] ..." << endl;
    cout << endl;
    cout << "  Compute eigenvalues and eigenvectors of a symmetric sparse matrix." << endl;
    cout << "  The matrix is stored in MATRIX_FILE in Matrix-Market format." << endl;
    cout << "  This driver invokes ARPACK for a Lanczos algorithm w/ implicit restarting." << endl;
    cout << endl;
    cout << "  -help               display this help and exit" << endl;
    cout << "  -srand              set random seed as time(NULL)" << endl;
    cout << "  -srand=NUMBER       set random seed as NUMBER" << endl;
    cout << "                      If an initial Lanczos vector is is not provided," << endl;
    cout << "                      a random initial vector is generated." << endl;
    cout << "  -printeig=NUMBER    0 for not to print eigenvalues and eigenvectors" << endl;
    cout << "                      1 for to print eigenvalues only (default)" << endl;
    cout << "                      2 for to print both eigenvalues and eigenvectors" << endl;
    cout << "  -nev=NUMBER         number of eigenvalues desired" << endl;
    cout << "  -ncv=NUMBER         (maximum) number of Lanczos basis vectors" << endl;
    cout << "  -part=STRING        LA, largest algebraic, for largest eigenvalues" << endl;
    cout << "                      SA, smallest algebraic, for smallest eigenvalues" << endl;
    cout << "                      BE, both ends, one more from high end if nev is odd" << endl;
    cout << "                      LM, largest magnitude, for largest eigenvalues in" << endl;
    cout << "                          magnitude (default)" << endl;
    cout << "                      SM, smallest magnitude, for smallest eigenvalues in" << endl;
    cout << "                          magnitude" << endl;
    cout << "  -tol=VALUE          tolerance for convergence test (default is 0.0)" << endl;
    cout << "                      If tol<=0.0, it will be set by ARPACK the default" << endl;
    cout << "                      tolerance, which is the machine epsilon." << endl;
    cout << endl;
    exit(1);
}
