//============================================================================
// Estimation of number of eigenvalues of a symmetric matrix by a randomized
// algorithm based on polynomial filtering.
//
// Reference:
// "Filtered Conjugate Residual-type Algorithms with Applications",
// Y. Saad, SIAM J. Matrix Anal. and Appl., Vol. 28, pp. 845--870, 2006.
//============================================================================

#include <stdlib.h>  // for exit, atof, atoi, atol, srand, rand, etc.
#include <string.h>  // for memcpy, strcmp, strncmp, etc.
#include <time.h>    // for time_t, time, clock_t, clock, CLOCKS_PER_SEC, etc.
#include <iostream>  // for cout, cerr, endl, etc. (under namespace std)

#include "matkit.h"
#include "symtrieig.h"
#include "laneig.h"
#include "polyfilt.h"
#include "filtlan.h"

using std::cout;
using std::cerr;
using std::endl;

void printUsageAndExit(const char *cmd);


int main(int argc, char *argv[]) {
    char *cmd = argv[0];
    if (argc <= 2                                // argc<=2 means no information about number of random vectors is provided
    || !strcmp(argv[1], "-?") || !strcmp(argv[1], "-h") || !strcmp(argv[1], "-help") || !strcmp(argv[1], "--help"))
        printUsageAndExit(cmd);                  // print usage and exit

    // declare variables for input parameters
    mkIndex baseDeg = 10, polyDeg = 0;           // degrees of base filter and polynomial filter
    Real bound0 = -Basics::infinity;             // the desired eigenvalues are in the interval [bound0, bound1]
    Real bound1 = Basics::infinity;
    Real eigLowerBound = -Basics::infinity;      // lower bound of eigenvalues
    Real eigUpperBound = Basics::infinity;       // upper bound of eigenvalues

    // Lanczos options; in case the min or max eigenvalues are not provided, the Lanczos algorithm will be used for computing them
    mkIndex numIterForEigenRange = 20;           // number of Lanczos iterations to determine the eigen-range
                                                 // if the eigen-range is provided, then numIterForEigenRange is ineffective,
                                                 // since the Lanczos algorithm will not be used

    IntervalOptions intervalOpts;                // options to determine the (transition) intervals and therefore the polynomial filter

    // parse for mintrials, stride, and maxtrials
    // number of estimated number of eigenvalues is reported for trials = mintrials:stride:maxtrials
    char *sptr = argv[2];
    int mintrials = atoi(sptr);
    sptr = strchr(sptr, ':') + 1;
    if (sptr-argv[2] <= 0) {
        cerr << "Error: the second parameter for numbers of trials is not in the right format." << endl;
        cerr << "It should be MINTRIALS:STRIDE:MAXTRIALS, for example, 10:10:100." << endl;
        exit(1);
    }
    int stride = atoi(sptr);
    sptr = strchr(sptr, ':') + 1;
    int maxtrials = atoi(sptr);
    if (sptr-argv[2] <= 0) {
        cerr << "Error: the second parameter for numbers of trials is not in the right format." << endl;
        cerr << "It should be MINTRIALS:STRIDE:MAXTRIALS, for example, 10:10:100." << endl;
        exit(1);
    }

    // adjust some parameters for this particualr application, estimation of number of eigenvalues
    intervalOpts.maxOuterIter = 0;  // don't check the balance, not really important in this application
    intervalOpts.initialShiftStep = 0.05;
    intervalOpts.transitionIntervalRatio = 2.0*intervalOpts.initialShiftStep;  // for high-pass filters (and low-pass filters w/ conversion)
    intervalOpts.initialPlateau = 1.0-2.0*intervalOpts.initialShiftStep;       // for mid-pass filters

    // interval weights for the polynomial filter
    // the weight of j-th interval is intervalWeights(j)
    intervalOpts.intervalWeights.resize(5);
    intervalOpts.intervalWeights(1) = 100.0;
    intervalOpts.intervalWeights(2) = 1.0;
    intervalOpts.intervalWeights(3) = 200.0;
    intervalOpts.intervalWeights(4) = 1.0;
    intervalOpts.intervalWeights(5) = 100.0;

    // parse input parameters
    for (int i=3; i<argc; i++) {
        if (!strcmp(argv[i], "-srand"))
            srand(time(NULL));                   // set random seed as time(NULL) for the random vector
        else if (!strncmp(argv[i], "-srand=", 7))
            srand(atoi(argv[i]+7));              // random seed (default 1) for the random vector

        else if (!strncmp(argv[i], "-bound0=", 8))
            bound0 = atof(argv[i]+8);            // lower bound of desired eigenvalues
        else if (!strncmp(argv[i], "-bound1=", 8))
            bound1 = atof(argv[i]+8);            // upper bound of desired eigenvalues
        else if (!strncmp(argv[i], "-basedeg=", 9))
            baseDeg = atoi(argv[i]+9);           // base filter degree
        else if (!strncmp(argv[i], "-polydeg=", 9))
            polyDeg = atoi(argv[i]+9);           // polynomial filter degree

        // if any of the two following parameters are not assigned, the standard Lanczos algorithm will be used to determine them
        else if (!strncmp(argv[i], "-mineig=", 8))
            eigLowerBound = atof(argv[i]+8);     // the smallest eigenvalue
        else if (!strncmp(argv[i], "-maxeig=", 8))
            eigUpperBound = atof(argv[i]+8);     // the largest eigenvalue

        // the following parameter numIterForEigenRange is effective only if the eigen-range (i.e. min/max eigenvalues) are not provided,
        // in which case the Lanczos algorithm will be invoked to approximate the two extreme eigenvalues
        else if (!strncmp(argv[i], "-eigrangeiter=", 14))
            numIterForEigenRange = atoi(argv[i]+14);            // number of Lanczos iterations for estimating the eigen-range

        // the following parameters are used for determine the (transition) intervals
        else if (!strncmp(argv[i], "-w1=", 4))
            intervalOpts.intervalWeights(1) = atof(argv[i]+4);  // weight of first  interval (default 100.0)
        else if (!strncmp(argv[i], "-w2=", 4))
            intervalOpts.intervalWeights(2) = atof(argv[i]+4);  // weight of second interval (default 1.0)
        else if (!strncmp(argv[i], "-w3=", 4))
            intervalOpts.intervalWeights(3) = atof(argv[i]+4);  // weight of third  interval (default 200.0)
        else if (!strncmp(argv[i], "-w4=", 4))
            intervalOpts.intervalWeights(4) = atof(argv[i]+4);  // weight of fourth interval (default 1.0, effective only for mid-pass filter)
        else if (!strncmp(argv[i], "-w5=", 4))
            intervalOpts.intervalWeights(5) = atof(argv[i]+4);  // weight of fifth  interval (default 100.0, effective only for mid-pass filter)
        else if (!strncmp(argv[i], "-ngrids=", 8))
            intervalOpts.numGridPoints = atoi(argv[i]+8);       // number of grid points
        else if (!strncmp(argv[i], "-delta=", 7))
            intervalOpts.initialShiftStep = atof(argv[i]+7);    // initial shift step relative to the length of interval of desired eigenvalues (default 0.05)
        else if (!strncmp(argv[i], "-maxouteriter=", 14))
            intervalOpts.maxOuterIter = atoi(argv[i]+14);       // maximum number of outer iterations to determine the (transition) intervals (default 0)
        else if (!strncmp(argv[i], "-ybottom=", 9))
            intervalOpts.yBottomLine = atof(argv[i]+9);         // the value which p(x) should be greater than for x in the interval of desired eigenvalues (default 0.001)
        else if (!strncmp(argv[i], "-yripple=", 9))
            intervalOpts.yRippleLimit = atof(argv[i]+9);        // the limit of height of ripples (not in the interval of desired eigenvalues) relative to the bottom of polynomial values for x in the interval of desired eigenvalues
        // the following parameter is effective only for high-pass filters (and low-pass filters w/ conversion)
        else if (!strncmp(argv[i], "-trans=", 7))
            intervalOpts.transitionIntervalRatio = atof(argv[i]+7);  // the (relative) length of transition interval
        // the following parameters are effective only for mid-pass filters
        else if (!strncmp(argv[i], "-plateau=", 9))
            intervalOpts.initialPlateau = atof(argv[i]+9);      // initial plateau relative to the length of interval in which eigenvalues are sought (default 0.1)
        else if (!strncmp(argv[i], "-ytol=", 6))
            intervalOpts.yLimitTol = atof(argv[i]+6);           // a mid-pass filter p(x) should have p(a1)=p(b1), where [a1,b1] is the interval of desired eigenvalues
                                                                // yLimitTol is the tolerance allowed for |p(a1)-p(b1)|
        else if (!strncmp(argv[i], "-maxinneriter=", 14))
            intervalOpts.maxInnerIter = atoi(argv[i]+14);       // maximum number of inner iterations to determine the (transition) intervals (default 30)
        else {
            cerr << "Error: the parameter \"" << argv[i] << "\" is not recognized!" << endl;
            cerr << "For more information, use \"" << cmd << " --help\"." << endl;
            exit(1);
        }
    }

    // deal with exceptions
    if (polyDeg == 0) {
        cerr << "Number of eigenvalues estimation driver: polynomial filter degree is not set!" << endl;
        exit(1);
    }
    if (bound0==-Basics::infinity && bound1==Basics::infinity) {
        cerr << "Number of eigenvalues estimation driver: the range of eigenvalues is not set!" << endl;
        exit(1);
    }
    if (bound0 >= bound1) {
        cerr << "Number of eigenvalues estimation driver: the lower bound (" << bound0 << ") should be less than the upper bound (" << bound1 << ")" << endl;
        exit(1);
    }

    // open the matrix file and read the matrix
    char *fileName = argv[1];
    cout << "# input matrix file: " << fileName << endl;
    SparseMatrix A = SparseMatrix::mmread(fileName);
    Real n = A.Nrows();
    cout << "# n=" << n << ", nnz=" << A.Nnz() << endl;

    // preprocessing
    if (eigLowerBound == -Basics::infinity || eigUpperBound == Basics::infinity) {
        if (eigLowerBound == -Basics::infinity && eigUpperBound == Basics::infinity)
            cout << "# the minimum and maximum eigenvalues are not provided, and therefore computed by the Lanczos algorithm" << endl;
        else if (eigLowerBound == -Basics::infinity)
            cout << "# the minimum eigenvalue is not provided, and therefore the Lanczos algorithm is invoked" << endl;
        else
            cout << "# the maximum eigenvalue is not provided, and therefore the Lanczos algorithm is invoked" << endl;
        clock_t start = clock();
        LanczosOptions opts;                  // default Lanczos options
        opts.wantEigVec = true;               // in case the default is false
        opts.minIter = numIterForEigenRange;  // want exactly numIterForEigenRange Lanczos iterations
        opts.maxIter = numIterForEigenRange;  // want exactly numIterForEigenRange Lanczos iterations
        Vector eigVal;
        Matrix eigVec;
        const char eigPart[] = "BE";
        LanczosEigenSolver(eigVal, eigVec, A, 2, eigPart, opts);
        eigLowerBound = eigVal(1) - (A*eigVec.column(1)-eigVal(1)*eigVec.column(1)).norm2();
        eigUpperBound = eigVal(2) + (A*eigVec.column(2)-eigVal(2)*eigVec.column(2)).norm2();
        cout << "# CPU time used for determining the eigen-range: " << (Real)(clock()-start) / (Real)CLOCKS_PER_SEC << " seconds" << endl;
    }
    if (bound1 <= eigLowerBound) {
        cerr << "the upper bound (" << bound1 << ") of eigenvalues to be sought should be greater than the lower bound of eigenvalues (" << eigLowerBound << ")!" << endl;
        exit(1);
    }
    if (bound0 >= eigUpperBound) {
        cerr << "the lower bound (" << bound0 << ") of eigenvalues to be sought should be less than the upper bound of eigenvalues (" << eigUpperBound << ")!" << endl;
        exit(1);
    }
    bound0 = Basics::max(bound0, eigLowerBound);
    bound1 = Basics::min(bound1, eigUpperBound);
    Vector frame(4);
    frame(1) = eigLowerBound;
    frame(2) = bound0;
    frame(3) = bound1;
    frame(4) = eigUpperBound;

    // print some input information
    cout << "# base filter degree = " << baseDeg << ", polynomial degree = " << polyDeg << endl;
    cout << "# all eigenvalues in the range: [" << frame(1) << ", " << frame(4) << "]" << endl;
    cout << "# eigenvalues of interest in the range: [" << frame(2) << ", " << frame(3) << "]" << endl;

    // estimate the number of eigenvalues in [bound0, bound1]
    Vector intervals = PolynomialFilterInterface::setFilter(A, frame, polyDeg, baseDeg, intervalOpts);
    cout << "# the partition for the filter: [" << intervals(1);
    for (mkIndex i=2; i<=intervals.Length(); i++)
        cout << ", " << intervals(i);
    cout << "]" << endl;
    Real meig = 0.0;
    for (int i=1; i<=maxtrials; i++) {
        Vector v = Vector::random(n, -1.0, 1.0);  // a random vector with elements uniformly distributed in [-1,1]
        v /= v.norm2();
        Vector w = PolynomialFilterInterface::filteredSparseMatrixPolynomialVectorProduct(v);  // p(A)*v
        meig += Vector::innerProduct(v, w);
        if (i >= mintrials && (i-mintrials)%stride == 0)
            cout << i << " " << meig*((Real)n/(Real)i) << endl;
    }

    return 0;
}

void printUsageAndExit(const char *cmd) {
    IntervalOptions opts0;
    cout << "Usage: " << cmd << " MATRIX_FILE NUM1:STRIDE:NUM2 [OPTION] [OPTION] ..." << endl;
    cout << endl;
    cout << "  Estimate the number of eigenvalues of a symmetric sparse matrix in a given" << endl;
    cout << "  interval. The matrix is stored in MATRIX_FILE in Matrix-Market format." << endl;
    cout << "  In total NUM2 random vectors will be generated, and the results are" << endl;
    cout << "  reported with from NUM1 to NUM2 random vectors for every STRIDE random" << endl;
    cout << "  vectors." << endl;
    cout << endl;
    cout << "  -help               display this help and exit" << endl;
    cout << "  -srand              set random seed as time(NULL)" << endl;
    cout << "  -srand=NUMBER       set random seed as NUMBER" << endl;
    cout << "                      The algorithm requires random vectors for estimating" << endl;
    cout << "                      the number of eigenvalues in the requested interval." << endl;
    cout << "  -bound0=VALUE       lower bound of eigenvalues to be sought" << endl;
    cout << "                      (default is no bound)" << endl;
    cout << "  -bound1=VALUE       upper bound of eigenvalues to be sought" << endl;
    cout << "                      (default is no bound)" << endl;
    cout << "  -mineig=VALUE       a (tight) lower bound of all eigenvalues" << endl;
    cout << "                      If not provided, it will be computed by the" << endl;
    cout << "                      Lanczos algorithm." << endl;
    cout << "  -maxeig=VALUE       a (tight) upper bound of all eigenvalues" << endl;
    cout << "                      If not provided, it will be computed by the" << endl;
    cout << "                      Lanczos algorithm." << endl;
    cout << "  -eigrangeiter=NUMBER  number of Lanczos iterations to determine the eigen" << endl;
    cout << "                      range, an interval which (tightly) contains all the" << endl;
    cout << endl;
    cout << "  Options for determining the polynomial filter:" << endl;
    cout << "  -basedeg=NUMBER     base filter degree (default 10) " << endl;
    cout << "                      basedeg from 8 to 15 usually work fine." << endl;
    cout << "  -polydeg=NUMBER     polynomial filter degree (suggested 10)" << endl;
    cout << "                      Note that polydeg should be increased if the requested" << endl;
    cout << "                      eigenvalues are interior especially in a narrow" << endl;
    cout << "                      interval relative to the range of the spectrum." << endl;
    cout << "  -w1=VALUE           weight of the first interval" << endl;
    cout << "  ..." << endl;
    cout << "  -w5=VALUE           weight of the fifth interval" << endl;
    cout << "                      (default w1,w2,w3,w4,w5 are " << opts0.intervalWeights(1) << ",";
    cout << opts0.intervalWeights(2) << "," << opts0.intervalWeights(3) << ",";
    cout << opts0.intervalWeights(4) << "," << opts0.intervalWeights(5) << ")" << endl;
    cout << "                      For high pass and low pass filters, there are 3" << endl;
    cout << "                      intervals, including 1 transition interval, so" << endl;
    cout << "                      only w1,w2,w3 are effective. A mid pass filter has 5" << endl;
    cout << "                      intervals, including 2 transition ones." << endl;
    cout << "  -ngrids=NUMBER      number of grid points, used to approximate the upper" << endl;
    cout << "                      bound of polynomial p(x) for x not in the interval in" << endl;
    cout << "                      in which the eigenvalues are requested" << endl;
    cout << "                      The value of ngrids (default " << opts0.numGridPoints << ") will automatically" << endl;
    cout << "                      be increased if it is too small." << endl;
    cout << "  -delta=VALUE        length of (initial) shift relative to the length of" << endl;
    cout << "                      the interval in which the eigenvalues are requested." << endl;
    cout << "                      The value should be between 0 and 1 (default " << opts0.initialShiftStep << ")." << endl;
    cout << "  -maxouteriter=NUMBER  maximum number of outer iterations to determine the" << endl;
    cout << "                      transition interval(s) (default is " << opts0.maxOuterIter << ")" << endl;
    cout << "  -ybottom=VALUE      the value which p(x) should be greater than for x in" << endl;
    cout << "                      [INTV2,INTV3] (default is " << opts0.yBottomLine << ")" << endl;
    cout << "  -yripple=VALUE      the limit of relative height of ripples not in the" << endl;
    cout << "                      interval of desired eigenvalues (default is " << opts0.yRippleLimit << ")" << endl;
    cout << "  -trans=VALUE        the relative length of transition interval, effective" << endl;
    cout << "                      only for high-pass filters (or low-pass filters with" << endl;
    cout << "                      conversion) (default is " << opts0.transitionIntervalRatio << ")" << endl;
    cout << endl;
    cout << "  The following parameters are effective only for mid-pass filters:" << endl;
    cout << "  -plateau=VALUE      length of (initial) plateau relative to the length of" << endl;
    cout << "                      the interval of desired eigenvalues" << endl;
    cout << "                      The value must be between 0 and 1 (default is " << opts0.initialPlateau << ")." << endl;
    cout << "  -ytol=VALUE         the tolerance allowed for |p(a1)-p(b1)|, where [a1,b1]" << endl;
    cout << "                      is the interval in which the eigenvalues are requested" << endl;
    cout << "                      (default is " << opts0.yLimitTol << ")" << endl;
    cout << "  -maxinneriter=NUMBER  maximum number of inner iterations to determine the" << endl;
    cout << "                      transition interval(s) (default is " << opts0.maxInnerIter << ")" << endl;
    cout << endl;
    exit(1);
}
