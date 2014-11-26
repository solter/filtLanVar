//============================================================================
// The laneig driver for computing extreme eigenvalues of a symmetric matrix
// by a shift-and-invert Lanczos procedure.
//
// coded by H.-r. Fang, last update May, 2012.
//============================================================================

#include <stdlib.h>  // for exit, atof, atoi, atol, srand, rand, etc.
#include <string.h>  // for memcpy, strcmp, strncmp, etc.
#include <time.h>    // for time_t, time, clock_t, clock, CLOCKS_PER_SEC, etc.
#include <iostream>  // for cout, cerr, endl, etc. (under namespace std)
#include <iomanip>   // for setw, setprecision, etc. (under namespace std)

#include "matkit.h"
#include "laneig.h"
#include "silaneig.h"

using std::cerr;
using std::cout;
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
    if (argc <= 2                                // argc<=2 means no input matrix is specified
    || !strcmp(argv[1], "-?") || !strcmp(argv[1], "-h") || !strcmp(argv[1], "-help") || !strcmp(argv[1], "--help"))
        printUsageAndExit(cmd);                  // print usage and exit

    // declare variables for input parameters
    bool printEigval = true;                     // whether to print eigenvalues or not
    bool printEigvec = false;                    // whether to print eigenvectors or not
    mkIndex neigWanted = 0;                      // number of eigenvalues desired
                                                 // default 0 means to obtain all eigenvalues in the specified range
    Real eigShift = Basics::infinity;            // the eigenvalues to be shifted by eigShift, i.e., A-eigShift*I
    char defaultEigPart[] = "LA";                // default part of desired eigenvalues
    char *eigPart = defaultEigPart;              // a string specifying which part of eigenvalues is desired
    ShiftAndInvertLanczosOptions opts;           // Lanczos options; see laneig.h for more information

    // parse input parameters
    for (int i=2; i<argc; i++) {
        if (!strcmp(argv[i], "-srand"))
            srand(time(NULL));                   // set random seed as time(NULL) for the initial Lanczos vector
        else if (!strncmp(argv[i], "-srand=", 7))
            srand(atoi(argv[i]+7));              // random seed (default 1) for the initial Lanczos vector
        else if (!strncmp(argv[i], "-printeig=", 10)) {
            if (atoi(argv[i]+10))                // whether to print eigenvalue or not
                printEigval = true;
            else
                printEigval = false;
            if (atoi(argv[i]+10) >= 2)           // whether to print eigenvector or not
                printEigvec = true;
            else
                printEigvec = false;
        }
        else if (!strncmp(argv[i], "-nev=", 5))
            neigWanted = atoi(argv[i]+5);        // number of eigenvalues desired (default 6)
        else if (!strncmp(argv[i], "-shift=", 7))
            eigShift = atof(argv[i]+7);          // the eigenvalues to be shifted by eigShift, i.e., A-eigShift*I
                                                 // (sparse) LU factorization of A-eigShift*I will be computed once and only once by UMFPACK
                                                 // (A-eigShift*I)^(-1)*v is computed via the LU factorization in each Lanczos iteration
        else if (!strncmp(argv[i], "-part=", 6))
                    eigPart = argv[i]+6;         // this string specifies which part of eigenvalues are desired, all closest to eigShift
                                                 // "LA", larger algebraically (default)
                                                 // "SA", smaller algebraically
                                                 // "BE", both sides equally, one more from the larger side if nev is odd
        else if (!strncmp(argv[i], "-miniter=", 9))
            opts.minIter = atoi(argv[i]+9);      // minimum number of the Lanczos iterations
        else if (!strncmp(argv[i], "-maxiter=", 9))
            opts.maxIter = atoi(argv[i]+9);      // maximum number of the Lanczos iterations
        else if (!strncmp(argv[i], "-extraiter=", 11))
            opts.extraIter = atoi(argv[i]+11);   // extra number of the Lanczos iterations after convergence, to avoid missing eigenvalues
        else if (!strncmp(argv[i], "-stride=", 8))
            opts.stride = atoi(argv[i]+8);       // convergence is checked every "stride" number of Lanczos iterations
        else if (!strncmp(argv[i], "-tol=", 5))
            opts.tol = atof(argv[i]+5);          // tolerance for convergence test
        else if (!strncmp(argv[i], "-bound0=", 8))
            opts.eigLowerBound = atof(argv[i]+8);  // lower bound of desired eigenvalues, effective if the the part of eigenvalues requested is "SA" or "BE"
        else if (!strncmp(argv[i], "-bound1=", 8))
            opts.eigUpperBound = atof(argv[i]+8);  // upper bound of desired eigenvalues, effective if the the part of eigenvalues requested is "LA" or "BE"
        else if (!strncmp(argv[i], "-bound=", 7)) {
            opts.eigLowerBound = atof(argv[i]+7);  // set both bounds as the same value
            opts.eigUpperBound = opts.eigLowerBound;
        }
        else if (!strncmp(argv[i], "-disp=", 6))
            opts.disp = atoi(argv[i]+6);         // diagnostic information display level (0, 1, or 2)
        else if (!strncmp(argv[i], "-expand=", 8))
            opts.memoryExpansionFactor = atof(argv[i]+8);  // it specifies how large the memory should be expanded when the allocated memory
                                                           // is insufficient to store the Lanczos vectors (default is 1.2)
        else if (!strncmp(argv[i], "-reorth=", 8))
            opts.reorth = atoi(argv[i]+8);       // 0 for no reorthogonalization (reserved)
                                                 // 1 for partial reorthogonalization (default)
                                                 // 2 for full reorthogonalization
        else if (!strncmp(argv[i], "-gamma2=", 8))
            opts.doubleReorthGamma = atof(argv[i]+8);  // reorthogonalization is doubled in the iterations with nrm < doubleReorthGamma*nrm_old, where
                                                       // nrm_old and nrm are the norms of the latest Lanczos vectors before and after reorthogonalization
        else if (!strncmp(argv[i], "-gamma1=", 8))
            opts.localReorthGamma = atof(argv[i]+8);   // local reorthogonalization is performed if beta[j-1] < localReorthGamma*beta[j] or localReorthGamma >= 1.0
                                                       // where beta[j] is the latest beta and beta[j-1] is the second to the latest beta
        else if (!strncmp(argv[i], "-strategy=", 10))
            opts.partialReorthStrategy = atoi(argv[i]+10);  // effective only if reorth==1 (i.e. partial reorthogonalization)
                                                            // 0 for reorthogonalization against previous Lanczos vectors in groups; each group consists of
                                                            //   v(i),v(i+1),...,v(j) with omega all greater than eta, and there is v(k) with omega(k) < delta, i<=k<=j
                                                            // 1 for reorthogonalization against previous Lanczos vectors with omega < eta
                                                            // 2 for reorthogonalization against previous Lanczos vectors v(i),v(i+1),...,v(j) such that
                                                            //   i is the smallest index for omega(i) > eta with the smallest i, and
                                                            //   j is the largest index for omega(j) > eta
                                                            // 3 for reorthogonalization against all previous Lanczos vectors
                                                            // omega(i) is the estimated orthogonal error of v(i) against the previous Lanczos vectors v(1),...,v(i-1)
                                                            // by default, the value of eta is eps^(0.75) and the value of delta is eps^(0.5),
                                                            // where eps is the machine epsilon
        else if (!strcmp(argv[i], "-checkorth"))
            opts.checkReorth = true;             // to check whether semi-orthogonalization is preserved with partial reorthogonalization
                                                 // expensive, should not be used for large problems
        else {
            cerr << "Error: the parameter \"" << argv[i] << "\" is not recognized!" << endl;
            cerr << "For more information, use \"" << cmd << " --help\"." << endl;
            exit(1);
        }
    }

    // deal with exceptions
    if (eigShift == Basics::infinity) {
        cerr << "Shift-and-invert Lanczos driver: shift is not set!" << endl;
        exit(1);
    }
    if (neigWanted == 0) {
        if (strcmp(eigPart, "LA") || strcmp(eigPart, "la")) {
            if (opts.eigUpperBound == Basics::infinity) {
                cerr << "Shift-and-invert Lanczos driver: eigen part is \"LA\", but the number of eigenvalues and high eigen bound are not specified!" << endl;
                exit(1);
            }
        }
        else if (strcmp(eigPart, "SA") || strcmp(eigPart, "sa")) {
            if (opts.eigLowerBound == -Basics::infinity) {
                cerr << "Shift-and-invert Lanczos driver: eigen part is \"SA\", but the number of eigenvalues and low eigen bound are not specified!" << endl;
                exit(1);
            }
        }
        else if (strcmp(eigPart, "BE") || strcmp(eigPart, "be")) {
            if (opts.eigLowerBound == -Basics::infinity || opts.eigUpperBound == Basics::infinity) {
                cerr << "Shift-and-invert Lanczos driver: eigen part is \"BE\", but the number of eigenvalues and high/low eigen bound are not specified!" << endl;
                exit(1);
            }
        }
    }

    // open the matrix file and read the matrix
    char *fileName = argv[1];
    cout << "input matrix file: " << fileName << endl;
    SparseMatrix A = SparseMatrix::mmread(fileName);
    cout << "n=" << A.Nrows() << ", nnz=" << A.Nnz() << endl;

    // print some input information
    if (opts.reorth == 0)
        cout << "Shift-and-invert Lanczos w/o reorthogonalization:" << endl;
    else {
        if (opts.reorth == 1) {
            cout << "Shift-and-invert Lanczos w/ partial reorthogonalization (strategy " << opts.partialReorthStrategy << ")" << endl;
            if (opts.localReorthGamma <= 0.0)
                cout << "extended local reorthogonalization is not applied" << endl;
            else if (opts.localReorthGamma < 1.0)
                cout << "extended local reorthogonalization gamma: " << opts.localReorthGamma << endl;
            else
                cout << "extended local reorthogonalization is enforced at every iteration" << endl;
        }
        else if (opts.reorth == 2)
            cout << "Shift-and-invert Lanczos w/ full reorthogonalization" << endl;
        else {
            cerr << "Error: unknown reorthogonalization index: " << opts.reorth << endl;
            exit(1);
        }
        if (opts.doubleReorthGamma <= 0.0)
            cout << "double reorthogonalization is not applied" << endl;
        else if (opts.doubleReorthGamma < 1.0)
            cout << "double reorthogonalization gamma: " << opts.doubleReorthGamma << endl;
        else
            cout << "double reorthogonalization is enforced at every iteration" << endl;
    }
    if (neigWanted != 0)
        cout << "nev=" << neigWanted << ", ";
    cout << "shift=" << eigShift << ", eigPart=" << eigPart;

    if (opts.eigLowerBound > -Basics::infinity && (!strcmp(eigPart,"SA") || !strcmp(eigPart,"sa") || !strcmp(eigPart,"BE") || !strcmp(eigPart,"BE"))) {
        cout << ", eigLowerBound=" << opts.eigLowerBound;
    }
    if (opts.eigUpperBound < Basics::infinity && (!strcmp(eigPart,"LA") || !strcmp(eigPart,"la") || !strcmp(eigPart,"BE") || !strcmp(eigPart,"BE"))) {
        cout << ", eigUpperBound=" << opts.eigUpperBound;
    }
    cout << endl;

    // computing eigenvalues by the shift-and-invert Lanczos algorithm
    clock_t start = clock();
    time_t start2 = time(NULL);
    Vector lambda;
    Matrix V;
    ShiftAndInvertLanczosInfo outputInfo = ShiftAndInvertLanczosEigenSolver(lambda, V, A, eigShift, neigWanted, eigPart, opts);
    cout << "minIter=" << opts.minIter << ", maxIter=" << opts.maxIter << ", stride=" << opts.stride << ", extraIter=" << opts.extraIter << ", tol=" << opts.tol << endl;
    Real totalCpuTime = (Real)(clock()-start) / (Real)CLOCKS_PER_SEC;
    time_t wallClockTime = time(NULL) - start2;

    // report the result
    cout << "========================================================" << endl;
    cout << "iteration count: " << outputInfo.numIter << endl;
    cout << "reorthogonalization iteration count: " << outputInfo.reorthIterCount << endl;
    cout << "reorthogonalization vector count: " << outputInfo.reorthVectorCount << " (" << 100.0*outputInfo.reorthVectorRate << "\%)" << endl;
    if (outputInfo.doubleReorthIterCount > 0)
        cout << "double reorthogonalization iteration count: " << outputInfo.doubleReorthIterCount << endl;
    if (outputInfo.localReorthIterCount > 0)
        cout << "local reorthogonalization iteration count: " << outputInfo.localReorthIterCount << endl;
    unsigned long memoryForMatrixInBytes = A.Nnz()*(sizeof(Real)+sizeof(mkIndex)) + (A.Ncols()+1)*sizeof(mkIndex);
    cout << "memory for storing matrix: " << memoryForMatrixInBytes << " bytes (" << Real(memoryForMatrixInBytes)/(Real)(1024*1024) << " MB)" << endl;
    cout << "memory required by Lanczos: " << outputInfo.memoryForLanczosInBytes << " bytes (" << Real(outputInfo.memoryForLanczosInBytes)/(Real)(1024*1024) << " MB)" << endl;
    cout << "CPU time for LU factorization: " << outputInfo.LuCpuTime << " seconds" << endl;
    cout << "CPU time for (A-s*I)^{-1}*v: " << outputInfo.forNextKrylovVectorCpuTime << " seconds" << endl;
    cout << "CPU time used for reorthogonalization: " << outputInfo.reorthogonalizationCpuTime << " seconds" << endl;
    cout << "CPU time used for convergence check: " << outputInfo.convergenceCheckCpuTime << " seconds" << endl;
    cout << "CPU time used for obtaining eigenvectors: " << outputInfo.forEigenVectorsCpuTime << " seconds" << endl;
    cout << "total CPU time: " << totalCpuTime << " seconds" << endl;
    cout << "wall clock time: " << wallClockTime << " seconds" << endl;
    cout << lambda.Length() << " eigenvalues found." << endl;
    if (lambda.Length() > 0) {
        Matrix E = A*V-V*lambda.spdiag();
        cout << "||A*V-V*S||_F=" << E.normFrobenius() << endl;
        cout << "||V'*V-I||_F=" << (V.transpose()*V-SparseMatrix::eye(lambda.Length())).normFrobenius() << endl;
        #ifdef USE_SINGLE
            const unsigned nprec = 7, nwidth= 17;
        #else
            const unsigned nprec = 15, nwidth= 25;
        #endif
        if (printEigval) {
            cout << "========================================================" << endl;
            cout << scientific << setprecision(nprec);
            cout << "eig id          eigenvalue             residual norm" << endl;
            for (mkIndex i=1; i<=lambda.Length(); i++) {
                cout << setw(6) << i << setw(nwidth) << lambda(i) << setw(nwidth) << E.column(i).norm2() << endl;
            }
        }
        if (printEigvec) {
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
    }

    return 0;
}

void printUsageAndExit(char* cmd) {
    ShiftAndInvertLanczosOptions opts0;
    cout << "Usage: " << cmd << " MATRIX_FILE [OPTION] [OPTION] ..." << endl;
    cout << endl;
    cout << "  Compute eigenvalues and eigenvectors of a symmetric sparse matrix." << endl;
    cout << "  The matrix is stored in MATRIX_FILE in Matrix-Market format." << endl;
    cout << endl;
    cout << "  -help               display this help and exit" << endl;
    cout << "  -srand              set random seed as time(NULL)" << endl;
    cout << "  -srand=NUMBER       set random seed as NUMBER" << endl;
    cout << "                      If an initial Lanczos vector is is not provided," << endl;
    cout << "                      a random initial vector is generated." << endl;
    cout << "  -printeig=NUMBER    0 for not to print eigenvalues and eigenvectors" << endl;
    cout << "                      1 for to print eigenvalues only (default)" << endl;
    cout << "                      2 for to print both eigenvalues and eigenvectors" << endl;
    cout << "  -shift=VALUE        all the eigenvalues to be shifted by VALUE" << endl;
    cout << "                      A sparse LU factorization of A-VALUE*I is computed" << endl;
    cout << "                      once and only once by UMFPACK. Then (A-VALUE*I)^(-1)*v" << endl;
    cout << "                      is computed via the LU factorization in each Lanczos" << endl;
    cout << "                      iteration." << endl;
    cout << "  -nev=NUMBER         number of eigenvalues desired" << endl;
    cout << "                      Use \"-nev=0\" for all eigenvalues in the specified" << endl;
    cout << "                      range (default)." << endl;
    cout << "  -part=STRING        LA, larger algebraically than the shift (default)" << endl;
    cout << "                      SA, smaller algebraically than the shift" << endl;
    cout << "                      BE, both sides equally of the shift, one more from the" << endl;
    cout << "                          larger side if nev is odd" << endl;
    cout << "  -miniter=NUMBER     minimum number of Lanczos iterations" << endl;
    cout << "  -maxiter=NUMBER     maximum number of Lanczos iterations" << endl;
    cout << "  -extraiter=NUMBER   extra number of Lanczos iterations, to improve the" << endl;
    cout << "                      accuracy and to avoid missing eigenvalues" << endl;
    cout << "  -stride=NUMBER      convergence is checked every NUMBER Lanczos" << endl;
    cout << "                      iterations" << endl;
    cout << "  -tol=VALUE          tolerance for convergence test" << endl;
    cout << "                      (default is " << opts0.tol << ")" << endl;
    cout << "  -bound0=VALUE       lower bound of eigenvalues to be sought" << endl;
    cout << "                      (default is no bound)" << endl;
    cout << "  -bound1=VALUE       upper bound of eigenvalues to be sought" << endl;
    cout << "                      (default is no bound)" << endl;
    cout << "                      If nev is 0 and part is LA, then bound1 must be set." << endl;
    cout << "                      If nev is 0 and part is SA, then bound0 must be set." << endl;
    cout << "                      If nev is 0 and part is BE, then both bound0 and" << endl;
    cout << "                      bound1 must be set." << endl;
    cout << "  -disp=INDEX         diagnostic information display level (default is " << opts0.disp << ")" << endl;
    cout << "                      0 for no information, 1 for some information, and" << endl;
    cout << "                      2 for more information" << endl;
    cout << "  -expand=VALUE       specifies how large the memory should be expanded" << endl;
    cout << "                      when the allocated memory is insufficient to store" << endl;
    cout << "                      the Lanczos vectors (default is " << opts0.memoryExpansionFactor << ")" << endl;
    cout << "  -reorth=INDEX       0 for no reorthogonalization (reserved option)" << endl;
    cout << "                      1 for partial reorthogonalization";
    if (opts0.reorth == 1)
        cout << " (default)";
    cout << endl;
    cout << "                      2 for full reorthogonalization";
    if (opts0.reorth == 2)
        cout << " (default)";
    cout << endl;
    cout << "  -gamma2=VALUE       reorthogonalization is doubled in the iterations with" << endl;
    cout << "                      nrm < gamma2*nrm_old, where nrm_old and nrm are the" << endl;
    cout << "                      norms of the latest Lanczos vectors before and after" << endl;
    cout << "                      reorthogonalization (default gamma2 is " << opts0.doubleReorthGamma << ")" << endl;
    cout << "                      If gamma2<=0.0, reorthogonalization is never doubled." << endl;
    cout << "                      If gamma2>=1.0, reorthogonalization is always doubled." << endl;
    cout << endl;
    cout << "  Options for partial reorthogonalization:" << endl;
    cout << "  -gamma1=VALUE       local reorthogonalization is performed in the" << endl;
    cout << "                      iterations with beta[j-1] < gamma1*beta[j] or" << endl;
    cout << "                      gamma1 >= 1.0, where beta[j] is the latest beta and" << endl;
    cout << "                      beta[j-1] is the second to the latest beta" << endl;
    cout << "                      (default gamma1 is " << opts0.localReorthGamma << ")" << endl;
    cout << "                      If gamma1<=0.0, local reorthogonalization is disabled." << endl;
    cout << "                      If gamma1>=1.0, local reorthogonalization is always" << endl;
    cout << "                      performed." << endl;
    cout << "  -strategy=INDEX     0 for reorthogonalization against previous Lanczos" << endl;
    cout << "                        vectors in groups; each group consists of v(i)," << endl;
    cout << "                        v(i+1),...,v(j) with omega all greater than eta, and" << endl;
    cout << "                        there is v(k) with omega(k) > delta, i<=k<=j" << endl;
    cout << "                      1 for reorthogonalization against previous Lanczos" << endl;
    cout << "                        vectors v(i) with omega(i) > eta" << endl;
    cout << "                      2 for reorthogonalization against previous Lanczos" << endl;
    cout << "                        vectors v(i),v(i+1),...,v(j) such that i is the" << endl;
    cout << "                        smallest index for omega(i) < eta, and j is the" << endl;
    cout << "                        largest index for omega(j) < eta" << endl;
    cout << "                      3 for reorthogonalization against all previous Lanczos" << endl;
    cout << "                        vectors, the most expensive option" << endl;
    cout << "                      (default strategy is " << opts0.partialReorthStrategy << ")" << endl;
    cout << "                      Here omega(i) is the estimated orthogonal error of" << endl;
    cout << "                      v(i) against the previous Lanczos vectors. By default," << endl;
    cout << "                      the value of eta=eps^0.75 and delta=eps^0.5, where eps" << endl;
    cout << "                      is the machine epsilon." << endl;
    cout << "  -checkorth          to check whether semi-orthogonalization is preserved" << endl;
    cout << "                      with partial reorthogonalization" << endl;
    cout << "                      This option should not be set for large problems," << endl;
    cout << "                      since it is expensive (default is ";
    if (opts0.checkReorth)
        cout << "to check)." << endl;
    else
        cout << "not to check)." << endl;
    cout << endl;
    exit(1);
}
