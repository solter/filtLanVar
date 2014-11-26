//============================================================================
// The filtlan driver for computing extreme or interior eigenvalues of a
// symmetric matrix by a filtered Lanczos procedure.
//
// Reference:
// "A Filtered Lanczos Procedure for Extreme and Interior Eigenvalue Problems"
// H.-r. Fang and Y. Saad, University of Minnesota Technical Report, 2011.
//============================================================================

#include <stdlib.h>  // for exit, atof, atoi, atol, srand, rand, etc.
#include <string.h>  // for memcpy, strcmp, strncmp, etc.
#include <time.h>    // for time_t, time, clock_t, clock, CLOCKS_PER_SEC, etc.
#include <iostream>  // for cout, cerr, endl, etc. (under namespace std)
#include <iomanip>   // for setw, setprecision, etc. (under namespace std)

#include "matkit.h"
#include "symtrieig.h"
#include "laneig.h"
#include "polyfilt.h"
#include "filtlan.h"

using std::cout;
using std::cerr;
using std::endl;

using std::setw;
using std::setprecision;
using std::scientific;

void printUsageAndExit(const char *cmd);


int main(int argc, char *argv[]) {
    char *cmd = argv[0];
    if (argc <= 2                                // argc<=2 means no input matrix is specified
    || !strcmp(argv[1], "-?") || !strcmp(argv[1], "-h") || !strcmp(argv[1], "-help") || !strcmp(argv[1], "--help"))
        printUsageAndExit(cmd);                  // print usage and exit

    // declare variables for input parameters
    bool printEigVal = true;                     // whether to print eigenvalues or not
    bool printEigVec = false;                    // whether to print eigenvectors or not

    FilteredLanczosOptions opts;                 // filtered Lanczos options; see filtlan.h and laneig.h for more information
    mkIndex baseDeg = 10, polyDeg = 0;           // degrees of base filter (default 10) and polynomial filter
    Vector frame(2);                             // the interval of the desired eigenvalues is [frame(1), frame(2)]
    frame(1) = -Basics::infinity;
    frame(2) = Basics::infinity;

    // parse input parameters
    for (int i=2; i<argc; i++) {
        if (!strcmp(argv[i], "-srand"))
            srand(time(NULL));                   // set random seed as time(NULL) for the initial Lanczos vector
        else if (!strncmp(argv[i], "-srand=", 7))
            srand(atoi(argv[i]+7));              // random seed (default 1) for the initial Lanczos vector
        else if (!strncmp(argv[i], "-printeig=", 10)) {
            if (atoi(argv[i]+10))                // whether to print eigenvalue or not
                printEigVal = true;
            else
                printEigVal = false;
            if (atoi(argv[i]+10) >= 2)           // whether to print eigenvector or not
                printEigVec = true;
            else
                printEigVec = false;
        }
        else if (!strncmp(argv[i], "-nev=", 5))
            opts.neigWanted = atoi(argv[i]+5);   // number of eigenvalues desired
        else if (!strncmp(argv[i], "-eigrangeiter=", 14))
            opts.numIterForEigenRange = atoi(argv[i]+14);  // number of Lanczos iterations to determine the eigen-range
        else if (!strncmp(argv[i], "-miniter=", 9))
            opts.minIter = atoi(argv[i]+9);      // minimum number of Lanczos iterations
        else if (!strncmp(argv[i], "-maxiter=", 9))
            opts.maxIter = atoi(argv[i]+9);      // maximum number of Lanczos iterations
        else if (!strncmp(argv[i], "-extraiter=", 11))
            opts.extraIter = atoi(argv[i]+11);   // extra number of Lanczos iterations after convergence, to avoid missing eigenvalues
        else if (!strncmp(argv[i], "-stride=", 8))
            opts.stride = atoi(argv[i]+8);       // convergence is checked every "stride" Lanczos iterations
        else if (!strncmp(argv[i], "-tol=", 5))
            opts.tol = atof(argv[i]+5);          // tolerance for convergence test
        else if (!strncmp(argv[i], "-bound0=", 8))
            frame(1) = atof(argv[i]+8);          // lower bound of desired eigenvalues
        else if (!strncmp(argv[i], "-bound1=", 8))
            frame(2) = atof(argv[i]+8);          // upper bound of desired eigenvalues
        else if (!strncmp(argv[i], "-basedeg=", 9))
            baseDeg = atoi(argv[i]+9);           // base filter degree
        else if (!strncmp(argv[i], "-polydeg=", 9))
            polyDeg = atoi(argv[i]+9);           // polynomial filter degree
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
        // the following parameters are used for determine the (transition) intervals
        else if (!strncmp(argv[i], "-w1=", 4))
            opts.intervalOpts.intervalWeights(1) = atof(argv[i]+4);  // weight of first  interval (default 100.0)
        else if (!strncmp(argv[i], "-w2=", 4))
            opts.intervalOpts.intervalWeights(2) = atof(argv[i]+4);  // weight of second interval (default 1.0)
        else if (!strncmp(argv[i], "-w3=", 4))
            opts.intervalOpts.intervalWeights(3) = atof(argv[i]+4);  // weight of third  interval (default 1.0)
        else if (!strncmp(argv[i], "-w4=", 4))
            opts.intervalOpts.intervalWeights(4) = atof(argv[i]+4);  // weight of fourth interval (default 1.0, effective only for mid-pass filter)
        else if (!strncmp(argv[i], "-w5=", 4))
            opts.intervalOpts.intervalWeights(5) = atof(argv[i]+4);  // weight of fifth  interval (default 100.0, effective only for mid-pass filter)
        else if (!strncmp(argv[i], "-ngrids=", 8))
            opts.intervalOpts.numGridPoints = atoi(argv[i]+8);       // number of grid points
        else if (!strncmp(argv[i], "-delta=", 7))
            opts.intervalOpts.initialShiftStep = atof(argv[i]+7);    // initial shift step relative to the length of interval of desired eigenvalues (default 0.01)
        else if (!strncmp(argv[i], "-maxouteriter=", 14))
            opts.intervalOpts.maxOuterIter = atoi(argv[i]+14);       // maximum number of outer iterations to determine the (transition) intervals (default 50)
        else if (!strncmp(argv[i], "-ybottom=", 9))
            opts.intervalOpts.yBottomLine = atof(argv[i]+9);         // the value which p(x) should be greater than for x in the interval of desired eigenvalues (default 0.001)
        else if (!strncmp(argv[i], "-yripple=", 9))
            opts.intervalOpts.yRippleLimit = atof(argv[i]+9);        // the limit of height of ripples (not in the interval of desired eigenvalues) relative to the bottom of polynomial values for x in the interval of desired eigenvalues
        // the following parameter is effective only for high-pass filters (and low-pass filters w/ conversion)
        else if (!strncmp(argv[i], "-trans=", 7))
            opts.intervalOpts.transitionIntervalRatio = atof(argv[i]+7);  // the (relative) length of transition interval (default 0.6)
        // the following parameters are effective only for mid-pass filters
        else if (!strncmp(argv[i], "-plateau=", 9))
            opts.intervalOpts.initialPlateau = atof(argv[i]+9);      // initial plateau relative to the length of interval of desired eigenvalues (default 0.1)
        else if (!strncmp(argv[i], "-ytol=", 6))
            opts.intervalOpts.yLimitTol = atof(argv[i]+6);           // a mid-pass filter p(x) should have p(a1)=p(b1), where [a1,b1] is the interval of desired eigenvalues
                                                                     // yLimitTol is the tolerance allowed for |p(a1)-p(b1)|
        else if (!strncmp(argv[i], "-maxinneriter=", 14))
            opts.intervalOpts.maxInnerIter = atoi(argv[i]+14);       // maximum number of inner iterations to determine the (transition) intervals (default 30)
        else {
            cerr << "Error: the parameter \"" << argv[i] << "\" is not recognized!" << endl;
            cerr << "For more information, use \"" << cmd << " --help\"." << endl;
            exit(1);
        }
    }

    // deal with exceptions
    if (polyDeg == 0) {
        cout << "Filtered Lanczos driver: the polynomial filter degree is not set!" << endl;
        exit(1);
    }
    if (baseDeg == 0)
        baseDeg = 10;
    if (frame(1)==-Basics::infinity && frame(2)==Basics::infinity) {
        cout << "Filtered Lanczos driver: the range of desired eigenvalues is not set!" << endl;
        exit(1);
    }
    if (frame(1) >= frame(2)) {
        cout << "Number of eigenvalues estimation driver: the lower bound (" << frame(1) << ") should be less than the upper bound (" << frame(2) << ")" << endl;
        exit(1);
    }

    // open the matrix file and read the matrix
    char *fileName = argv[1];
    cout << "input matrix file: " << fileName << endl;
    SparseMatrix A = SparseMatrix::mmread(fileName);
    cout << "n=" << A.Nrows() << ", nnz=" << A.Nnz() << endl;

    // print some input information
    if (opts.reorth == 0)
        cout << "Filtered Lanczos w/o reorthogonalization" << endl;
    else if (opts.reorth == 1) {
       if (opts.reorth == 1) {
            cout << "Filtered Lanczos w/ partial reorthogonalization (strategy " << opts.partialReorthStrategy << ")" << endl;
            if (opts.localReorthGamma <= 0.0)
                cout << "extended local reorthogonalization is not applied" << endl;
            else if (opts.localReorthGamma < 1.0)
                cout << "extended local reorthogonalization gamma: " << opts.localReorthGamma << endl;
            else
                cout << "extended local reorthogonalization is enforced at every iteration" << endl;
        }
        else if (opts.reorth == 2)
            cout << "Filtered Lanczos w/ full reorthogonalization" << endl;
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
    if (opts.neigWanted != 0)  // default is 0, which means to obtain all eigenvalues in the specified interval
        cout << "number of eigenvalues / eigenvectors wanted: " << opts.neigWanted << endl;
    cout << "base filter degree: " << baseDeg << endl;
    cout << "polynomial degree: " << polyDeg << endl;
    cout << "eigenvalues desired in the range: [" << frame(1) << ", " << frame(2) << "]" << endl;

    // compute eigenvalues by the filtered Lanczos algorithm
    clock_t start = clock();
    time_t start2 = time(NULL);
    Vector lambda;
    Matrix V;
    FilteredLanczosInfo outputInfo = FilteredLanczosEigenSolver(lambda, V, A, frame, polyDeg, baseDeg, opts);
    cout << "minIter=" << opts.minIter << ", maxIter=" << opts.maxIter << ", stride=" << opts.stride << ", extraIter=" << opts.extraIter << ", tol=" << opts.tol << endl;
    Real totalCpuTime = (Real)(clock()-start) / (Real)CLOCKS_PER_SEC;
    time_t wallClockTime = time(NULL) - start2;

    // report the result
    cout << "========================================================" << endl;
    cout << "eigenvalues all in the range: [" << frame(1) << ", " << frame(4) << "]" << endl;
    cout << "interval weights: " << opts.intervalOpts.intervalWeights(1) << ", " <<  opts.intervalOpts.intervalWeights(2) << ", " << opts.intervalOpts.intervalWeights(3);
    if (frame(1) != frame(2) && frame(3) != frame(4))
        cout << ", " << opts.intervalOpts.intervalWeights(4) << ", " << opts.intervalOpts.intervalWeights(5);
    cout << endl;
    cout << "intervals: [" << outputInfo.intervals(1);
    for (mkIndex i=2; i<=outputInfo.intervals.Length(); i++)
        cout << ", " << outputInfo.intervals(i);
    cout << "]" << endl;
    cout << "iteration count: " << outputInfo.numIter << endl;
    cout << "reorthogonalization iteration count: " << outputInfo.reorthIterCount << endl;
    cout << "reorthogonalization vector count: " << outputInfo.reorthVectorCount << " (" << 100.0*outputInfo.reorthVectorRate << "%)" << endl;
    if (outputInfo.doubleReorthIterCount > 0)
        cout << "double reorthogonalization iteration count: " << outputInfo.doubleReorthIterCount << endl;
    if (outputInfo.localReorthIterCount > 0)
        cout << "extended local reorthogonalization iteration count: " << outputInfo.localReorthIterCount << endl;
    unsigned long memoryForMatrixInBytes = A.Nnz()*(sizeof(Real)+sizeof(mkIndex)) + (A.Ncols()+1)*sizeof(mkIndex);
    cout << "memory for storing matrix: " << memoryForMatrixInBytes << " bytes (" << Real(memoryForMatrixInBytes)/(Real)(1024*1024) << " MB)" << endl;
    cout << "memory required by Lanczos: " << outputInfo.memoryForLanczosInBytes << " bytes (" << Real(outputInfo.memoryForLanczosInBytes)/(Real)(1024*1024) << " MB)" << endl;
    cout << "CPU time used for determining the eigen-range: " << outputInfo.forEigenRangeCpuTime << " seconds" << endl;
    cout << "(# iter: " << outputInfo.forEigenRangeInfo.numIter << ";";
    cout << " matvec time: " << outputInfo.forEigenRangeInfo.forNextKrylovVectorCpuTime << " seconds;";
    cout << " reorth time: " << outputInfo.forEigenRangeInfo.reorthogonalizationCpuTime << " seconds;";
    cout << " convergence check time: " << outputInfo.forEigenRangeInfo.convergenceCheckCpuTime << " seconds)"<< endl;
    cout << "CPU time for p(A)*v: " << outputInfo.forNextKrylovVectorCpuTime << " seconds" << endl;
    cout << "CPU time used for reorthogonalization: " << outputInfo.reorthogonalizationCpuTime << " seconds" << endl;
    cout << "CPU time used for convergence check: " << outputInfo.convergenceCheckCpuTime << " seconds" << endl;
    cout << "CPU time used for obtaining eigenvectors: " << outputInfo.forEigenVectorsCpuTime << " seconds" << endl;
    cout << "total CPU time used: " << totalCpuTime << " seconds" << endl;
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
    }

    return 0;
}

void printUsageAndExit(const char *cmd) {
    FilteredLanczosOptions opts0;
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
    cout << "  -nev=NUMBER         number of eigenvalues desired" << endl;
    cout << "  -eigrangeiter=NUMBER  number of Lanczos iterations to determine the eigen" << endl;
    cout << "                      range, an interval which (tightly) contains all the" << endl;
    cout << "                      eigenvalues (default is " << opts0.numIterForEigenRange << ")" << endl;
    cout << "  -miniter=NUMBER     minimum number of Lanczos iterations" << endl;
    cout << "  -maxiter=NUMBER     maximum number of Lanczos iterations" << endl;
    cout << "  -extraiter=NUMBER   extra number of Lanczos iterations, to avoid missing" << endl;
    cout << "                      eigenvalues and to improve the accuracy" << endl;
    cout << "                      (default is " << opts0.extraIter << ")" << endl;
    cout << "  -stride=NUMBER      convergence is checked every NUMBER Lanczos" << endl;
    cout << "                      iterations" << endl;
    cout << "  -tol=VALUE          tolerance for convergence test" << endl;
    cout << "                      (default is " << opts0.tol << ")" << endl;
    cout << "  -bound0=VALUE       lower bound of eigenvalues to be sought" << endl;
    cout << "                      (default is no bound)" << endl;
    cout << "  -bound1=VALUE       upper bound of eigenvalues to be sought" << endl;
    cout << "                      (default is no bound)" << endl;
    cout << "                      At least one of bound0 and bound1 must be set. If only" << endl;
    cout << "                      one bound is set, it is an extreme eigenvalue problem." << endl;
    cout << "                      If both bounds are set, it is an interior eigenvalue" << endl;
    cout << "                      problem." << endl;
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
    cout << "  Options for determining the polynomial filter:" << endl;
    cout << "  -basedeg=NUMBER     base filter degree (suggested 10) " << endl;
    cout << "                      basedeg from 8 to 15 usually work fine." << endl;
    cout << "  -polydeg=NUMBER     polynomial filter degree (suggested 10)" << endl;
    cout << "                      Note that polydeg should be increased if the requested" << endl;
    cout << "                      eigenvalues are interior especially in a narrow" << endl;
    cout << "                      interval relative to the range of the spectrum." << endl;
    cout << "  -w1=VALUE           weight of the first interval" << endl;
    cout << "  ..." << endl;
    cout << "  -w5=VALUE           weight of the fifth interval" << endl;
    cout << "                      (default w1,w2,w3,w4,w5 are " << opts0.intervalOpts.intervalWeights(1) << ",";
    cout << opts0.intervalOpts.intervalWeights(2) << "," << opts0.intervalOpts.intervalWeights(3) << ",";
    cout << opts0.intervalOpts.intervalWeights(4) << "," << opts0.intervalOpts.intervalWeights(5) << ")" << endl;
    cout << "                      For high pass and low pass filters, there are 3" << endl;
    cout << "                      intervals, including 1 transition interval, so" << endl;
    cout << "                      only w1,w2,w3 are effective. A mid pass filter has 5" << endl;
    cout << "                      intervals, including 2 transition ones." << endl;
    cout << "  -ngrids=NUMBER      number of grid points, used to approximate the upper" << endl;
    cout << "                      bound of polynomial p(x) for x not in the interval in" << endl;
    cout << "                      which the eigenvalues are sought" << endl;
    cout << "                      The value of ngrids (default " << opts0.intervalOpts.numGridPoints << ") will automatically" << endl;
    cout << "                      be increased if it is too small." << endl;
    cout << "  -delta=VALUE        length of (initial) shift relative to the length of" << endl;
    cout << "                      the interval in which the eigenvalues are requested." << endl;
    cout << "                      The value should be between 0 and 1 (default " << opts0.intervalOpts.initialShiftStep << ")." << endl;
    cout << "  -maxouteriter=NUMBER  maximum number of outer iterations to determine the" << endl;
    cout << "                      transition interval(s) (default is " << opts0.intervalOpts.maxOuterIter << ")" << endl;
    cout << "  -ybottom=VALUE      the value which p(x) should be greater than for x in" << endl;
    cout << "                      [INTV2,INTV3] (default is " << opts0.intervalOpts.yBottomLine << ")" << endl;
    cout << "  -yripple=VALUE      the limit of relative height of ripples not in the" << endl;
    cout << "                      interval of desired eigenvalues (default is " << opts0.intervalOpts.yRippleLimit << ")" << endl;
    cout << "  -trans=VALUE        the relative length of transition interval, effective" << endl;
    cout << "                      only for high-pass filters (or low-pass filters with" << endl;
    cout << "                      conversion) (default is " << opts0.intervalOpts.transitionIntervalRatio << ")" << endl;
    cout << endl;
    cout << "  The following parameters are effective only for mid-pass filters:" << endl;
    cout << "  -plateau=VALUE      length of (initial) plateau relative to the length of" << endl;
    cout << "                      the interval of desired eigenvalues" << endl;
    cout << "                      The value must be between 0 and 1 (default is " << opts0.intervalOpts.initialPlateau << ")." << endl;
    cout << "  -ytol=VALUE         the tolerance allowed for |p(a1)-p(b1)|, where [a1,b1]" << endl;
    cout << "                      is the interval in which the eigenvalues are requested" << endl;
    cout << "                      (default is " << opts0.intervalOpts.yLimitTol << ")" << endl;
    cout << "  -maxinneriter=NUMBER  maximum number of inner iterations to determine the" << endl;
    cout << "                      transition interval(s) (default is " << opts0.intervalOpts.maxInnerIter << ")" << endl;
    cout << endl;
    exit(1);
}
