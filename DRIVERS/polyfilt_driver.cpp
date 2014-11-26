//============================================================================
// A driver for generating plotting data of a polynomial filter
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
#include "polyfilt.h"

using std::cout;
using std::cerr;
using std::endl;

void printUsageAndExit(const char *cmd);


int main(int argc, char *argv[]) {
    char *cmd = argv[0];
    if (argc < 4                                // means no sufficient arguments
        || !strcmp(argv[1], "-?") || !strcmp(argv[1], "-h") || !strcmp(argv[1], "-help") || !strcmp(argv[1], "--help"))
        printUsageAndExit(cmd);                 // print usage and exit

    // base filter degree and polynomial degree
    mkIndex baseDeg = atoi(argv[1]);
    mkIndex polyDeg = atoi(argv[2]);

    // initial intervals
    char *sptr = argv[3];
    Vector frame(4);
    frame(1) = atof(sptr);
    sptr = strchr(sptr, ',') + 1;
    if (sptr-argv[3] <= 0) {
        cerr << "Error: the third parameter, the intervals, is not in the right format." << endl;
        cerr << "It should be INTV1,INTV2,INTV3,INTV4 such that INTV1<=INTV2<INTV3<=INTV4;" << endl;
        cerr << "all the eigenvalues are (tightly) in [INTV1,INTV4], and the desired" << endl;
        cerr << "eigenvalues are in [INTV2,INTV3]." << endl;
        exit(1);
    }
    frame(2) = atof(sptr);
    sptr = strchr(sptr, ',') + 1;
    if (sptr-argv[3] <= 0) {
        cerr << "Error: the third parameter, the intervals, is not in the right format." << endl;
        cerr << "It should be INTV1,INTV2,INTV3,INTV4 such that INTV1<=INTV2<INTV3<=INTV4;" << endl;
        cerr << "all the eigenvalues are (tightly) in [INTV1,INTV4], and the desired" << endl;
        cerr << "eigenvalues are in [INTV2,INTV3]." << endl;
        exit(1);
    }
    frame(3) = atof(sptr);
    sptr = strchr(sptr, ',') + 1;
    if (sptr-argv[3] <= 0) {
        cerr << "Error: the third parameter, the intervals, is not in the right format." << endl;
        cerr << "It should be INTV1,INTV2,INTV3,INTV4 such that INTV1<=INTV2<INTV3<=INTV4;" << endl;
        cerr << "all the eigenvalues are (tightly) in [INTV1,INTV4], and the desired" << endl;
        cerr << "eigenvalues are in [INTV2,INTV3]." << endl;
        exit(1);
    }
    frame(4) = atof(sptr);

    // intervals are refined by adding transition interval(s), up to 6 intervals / interval weights
    Vector intervalWeights = Vector::ones(5);
    intervalWeights(1) = 200.0;
    intervalWeights(5) = 200.0;

    // interval options
    IntervalOptions opts;

    // number of points for graph plot
    mkIndex npoints = 1000;
    npoints = Basics::max(npoints, (mkIndex)(8.0*(frame(4)-frame(1))/(frame(3)-frame(2))));

    // parse input parameters
    for (int i=4; i<argc; i++) {
        if (!strncmp(argv[i], "-npoints=", 9))
            npoints = atoi(argv[i]+9);                // number of graph points
        else if (!strncmp(argv[i], "-w1=", 4))
            intervalWeights(1) = atof(argv[i]+4);     // weight of first  interval (default 200.0)
        else if (!strncmp(argv[i], "-w2=", 4))
            intervalWeights(2) = atof(argv[i]+4);     // weight of second interval (default 1.0)
        else if (!strncmp(argv[i], "-w3=", 4))
            intervalWeights(3) = atof(argv[i]+4);     // weight of third  interval (default 1.0)
        else if (!strncmp(argv[i], "-w4=", 4))
            intervalWeights(4) = atof(argv[i]+4);     // weight of fourth interval (default 1.0, effective only for mid-pass filter)
        else if (!strncmp(argv[i], "-w5=", 4))
            intervalWeights(5) = atof(argv[i]+4);     // weight of fifth  interval (default 200.0, effective only for mid-pass filter)
        else if (!strncmp(argv[i], "-ngrids=", 8))
            opts.numGridPoints = atoi(argv[i]+8);     // number of grid points
        else if (!strncmp(argv[i], "-delta=", 7))
            opts.initialShiftStep = atof(argv[i]+7);  // initial shift step relative to the length of interval of desired eigenvalues (default 0.01)
        else if (!strncmp(argv[i], "-maxouteriter=", 14))
            opts.maxOuterIter = atoi(argv[i]+14);     // maximum number of outer iterations to determine the (transition) intervals (default 50)
        else if (!strncmp(argv[i], "-ybottom=", 9))
            opts.yBottomLine = atof(argv[i]+9);       // the value which p(x) should be greater than for x in the interval of desired eigenvalues (default 0.001)
        else if (!strncmp(argv[i], "-yripple=", 9))
            opts.yRippleLimit = atof(argv[i]+9);      // the limit of height of ripples (not in the interval of desired eigenvalues) relative to the bottom of polynomial values for x in the interval of desired eigenvalues
        // the following parameter is effective only for high-pass filters (and low-pass filters w/ conversion)
        else if (!strncmp(argv[i], "-trans=", 7))
            opts.transitionIntervalRatio = atof(argv[i]+7);  // the (relative) length of transition interval (default 0.6)
        // the following parameters are effective only for mid-pass filters
        else if (!strcmp(argv[i], "-reverseintv"))
            opts.reverseInterval = true;              // whether to reverse the interval or not
        else if (!strncmp(argv[i], "-plateau=", 9))
            opts.initialPlateau = atof(argv[i]+9);    // initial plateau relative to the length of interval of desired eigenvalues (default 0.1)
        else if (!strncmp(argv[i], "-ytol=", 6))
            opts.yLimitTol = atof(argv[i]+6);         // a mid-pass filter p(x) should have p(a1)=p(b1), where [a1,b1] is the interval of desired eigenvalues
                                                      // yLimitTol is the tolerance allowed for |p(a1)-p(b1)|
        else if (!strncmp(argv[i], "-maxinneriter=", 14))
            opts.maxInnerIter = atoi(argv[i]+14);     // maximum number of inner iterations to determine the (transition) intervals (default 30)
        else {
            cerr << "Error: the parameter \"" << argv[i] << "\" is not recognized!" << endl;
            cerr << "For more information, use \"" << cmd << " --help\"." << endl;
            exit(1);
        }
    }

    // print some information
    cout << "# base filter degree: " << baseDeg << endl;
    cout << "# polynomial degree: " << polyDeg << endl;
    cout << "# frames: " << frame(1) << ", " << frame(2) << ", " << frame(3) << ", " << frame(4) << endl;
    cout << "# interval weights: " << intervalWeights(1) << ", " <<  intervalWeights(2) << ", " << intervalWeights(3);
    if (frame(1) != frame(2) && frame(3) != frame(4))
        cout << ", " << intervalWeights(4) << ", " << intervalWeights(5);
    cout << endl;
    if (frame(1) != frame(2) && frame(3) != frame(4)) {
        cout << "# maximum number of outer iterations: " << opts.maxOuterIter << endl;
        cout << "# maximum number of inner iterations: " << opts.maxInnerIter << endl;
    }
    else
        cout << "# maximum number of iterations: " << opts.maxOuterIter << endl;

    // reverse-and-shift the interval, if required
    bool intervalReversed = opts.reverseInterval;  // 1.0 means not to reverse; -1.0 means to reverse
    if (frame(1) == frame(2))
        intervalReversed = true;                   // a low-pass filter is expected for smallest eigenvalues; reverse the interval and find a high-pass filter instead
    else if (frame(3) == frame(4))
        intervalReversed = false;                  // interval cannot be reversed for a high-pass filter
    Vector frame2;
    if (intervalReversed)
        frame2 = -frame.reverse();
    else
        frame2 = frame;
    Real shift = frame2(1);
    frame2 -= shift;

    // get the intervals, i.e., adding transition interval(s)
    Vector intervals2;
    clock_t start = clock();
    PolynomialFilterInfo filterInfo = GetIntervals(intervals2, frame2, polyDeg, baseDeg, opts);
    Real cpuTime = (clock()-start) / (Real)CLOCKS_PER_SEC;
    Vector intervals = intervals2 + shift;
    if (intervalReversed)
        intervals = -intervals.reverse();

    // report the result
    if (filterInfo.filterOK == 2)
        cout << "# the filter is optimal (quality index: " << filterInfo.filterQualityIndex << ")" << endl;
    else if (filterInfo.filterOK == 1)
        cout << "# the filter is OK (quality index: " << filterInfo.filterQualityIndex << ")" << endl;
    else
        cout << "# the filter is no good" << endl;
    cout << "# CPU time to get the interval / filter: " << cpuTime << " seconds" << endl;
    cout << "# interval: " << intervals(1);
    for (mkIndex i=2; i<=intervals.Length(); i++)
        cout << ", " << intervals(i);
    cout << endl;
    cout << "# y-summit: " << filterInfo.ySummit << endl;
    cout << "# y-limit: " << filterInfo.yLimit << endl;
    if (frame2(3) == frame2(4)) {
        cout << "# y-summit (wing): " << filterInfo.yLeftSummit << endl;
        cout << "# y-bottom (wing): " << filterInfo.yLeftBottom << endl;
    }
    else {
        cout << "# y-limit-gap: " << filterInfo.yLimitGap << endl;
        if (intervalReversed) {
            cout << "# y-summit (left wing): " << filterInfo.yRightSummit << endl;
            cout << "# y-bottom (left wing): " << filterInfo.yRightBottom << endl;
            cout << "# y-summit (right wing): " << filterInfo.yLeftSummit << endl;
            cout << "# y-bottom (right wing): " << filterInfo.yLeftBottom << endl;
        }
        else {
            cout << "# y-summit (left wing): " << filterInfo.yLeftSummit << endl;
            cout << "# y-bottom (left wing): " << filterInfo.yLeftBottom << endl;
            cout << "# y-summit (right wing): " << filterInfo.yRightSummit << endl;
            cout << "# y-bottom (right wing): " << filterInfo.yRightBottom << endl;
        }
    }
    cout << "# number of grid points: " << opts.numGridPoints << endl;
    cout << "# no. iterations to get the filter: " << filterInfo.numIter << endl;
    cout << "# total no. iterations: " << filterInfo.totalNumIter << endl;
    if (frame2(3) != frame2(4)) {
        if (intervalReversed) {
            cout << "# no. left steps: " << filterInfo.numRightSteps << endl;
            cout << "# no. right steps: " << filterInfo.numLeftSteps << endl;
        }
        else {
            cout << "# no. left steps: " << filterInfo.numLeftSteps << endl;
            cout << "# no. right steps: " << filterInfo.numRightSteps << endl;
        }
    }
    cout << "# number of graph points: " << npoints << endl;

    // output the graph points
    const int HighLowFlags[] = { 1, -1, 0, -1, 1 };
    Matrix baseFilter = HermiteBaseFilterInChebyshevBasis(intervals2, HighLowFlags, baseDeg);
    Matrix polyFilter = FilteredConjugateResidualPolynomial(baseFilter, intervals2, intervalWeights, polyDeg);
    Vector yb(npoints), yp(npoints);
    Real step = (frame(4)-frame(1)) / (Real)(npoints-1);
    Real x = frame2(1);  // which is 0.0
    if (intervalReversed) {
        for (mkIndex i=0; i<npoints; i++) {
            yb(npoints-i) = 1.0 - PiecewisePolynomialEvaluationInChebyshevBasis(baseFilter, intervals2, x);
            yp(npoints-i) = 1.0 - PiecewisePolynomialEvaluationInChebyshevBasis(polyFilter, intervals2, x);
            x += step;
        }
    }
    else {
        for (mkIndex i=1; i<=npoints; i++) {
            yb(i) = 1.0 - PiecewisePolynomialEvaluationInChebyshevBasis(baseFilter, intervals2, x);
            yp(i) = 1.0 - PiecewisePolynomialEvaluationInChebyshevBasis(polyFilter, intervals2, x);
            x += step;
        }
    }
    x = frame(1);
    #ifdef USE_SINGLE
        const unsigned nprec = 7, nwidth= 17;
    #else
        const unsigned nprec = 15, nwidth= 25;
    #endif
    cout << std::scientific << std::setprecision(nprec);
    for (mkIndex i=1; i<=npoints; i++) {
        cout << std::setw(nwidth) << x << std::setw(nwidth) << filterInfo.yLimit << std::setw(nwidth) << yb(i) << std::setw(nwidth) << yp(i) << endl;
        x += step;
    }

    return 0;
}

void printUsageAndExit(const char *cmd) {
    IntervalOptions opts0;
    cout << "Usage: " << cmd << " BASEDEG POLYDEG IT1,IT2,IT3,IT4 [OPTION] [OPTION] ..." << endl;
    cout << endl;
    cout << "  Generate data for plotting a polynomial filter (e.g. by gnuplot)." << endl;
    cout << endl;
    cout << "  BASEDEG is the base filter degree." << endl;
    cout << "  POLYDEG is the polynomial degree." << endl;
    cout << endl;
    cout << "  [IT1, IT4] is a (tight) interval containing all eigenvalues." << endl;
    cout << "  [IT2, IT3] is the interval in which the eigenvalues are of interest." << endl;
    cout << "  IT1 <= IT2 < IT3 <= IT4" << endl;
    cout << "  If INTV1==INTV2, then INTV3<INTV4 and the polynomial is a low pass filter." << endl;
    cout << "  If INTV3==INTV4, then INTV1<INTV2 and the polynomial is a high pass filter." << endl;
    cout << "  Otherwise, the polynomial is a mid pass filter." << endl;
    cout << endl;
    cout << "  Options:" << endl;
    cout << "  -npoints=NUMBER     number of graphing points which determines the plot" << endl;
    cout << "                      resolution (default is 1000)" << endl;
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
