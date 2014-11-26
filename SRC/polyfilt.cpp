//============================================================================
// Routines for polynomial approximation, Chebyshev polynomial computation, a
// conjugate-residual-type algorithm in polynomial space, etc.
//
// Reference:
// "A Filtered Lanczos Procedure for Extreme and Interior Eigenvalue Problems"
// H.-r. Fang and Y. Saad, University of Minnesota Technical Report, 2011.
//============================================================================

#include <stdlib.h>  // for exit, atof, atoi, atol, srand, rand, etc.
#include <string.h>  // for memcpy, strcmp, strncmp, etc.
#include <stddef.h>  // for NULL pointer, pointer subtraction, etc.
#include <math.h>    // for sqrt, pow, log10, etc.
#include <iostream>  // for cout, cerr, endl, etc. (under namespace std)

#include "matkit.h"
#include "polyfilt.h"

using std::endl;

#define USE_CSRMV
// A*v corresponds to CSCMV and v*A corresponds to CSRMV
// when A is (sparse) symmetric, then A*v equals v*A in exact arithmetic
// however CSRMV usually slightly faster than CSCMV
// "#define USE_CSRMV" invokes CSRMV


////////////////////////////////////////////////////////////////////////////
//    Newton - Hermite Polynomial Interpolation
////////////////////////////////////////////////////////////////////////////

// build P(z) by Newton's divided differences in the form
//     P(z) = a(1) + a(2)*(z-x(1)) + a(3)*(z-x(1))*(z-x(2)) + ... + a(n)*(z-x(1))*...*(z-x(n-1)),
// such that P(x(i)) = y(i) for i=1,...,n, where
//     x,y are input vectors of length n, and a is the output vector of length n
// if x(i)==x(j) for some i!=j, then it is assumed that the derivative of P(z) is to be zero at x(i),
//     and the Hermite polynomial interpolation is applied
// in general, if there are k x(i)'s with the same value x0, then
//     the j-th order derivative of P(z) is zero at z=x0 for j=1,...,k-1
Vector NewtonPolynomial(const Vector &x, const Vector &y) {
    mkIndex n = x.Length();
    if (n != y.Length()) {
        *Basics::err << "NewtonPolynomial(const Vector &x, const Vector &y): x and y must have the same number of entries, i.e. x.Length()==y.Length()!" << endl;
        Basics::quit(1);
    }
    if (n == 0u) {
        // exception, return a zero vector
        return Vector();
    }

    // input storage
    const Real *sx = x.Store();
    const Real *sy = y.Store();

    // output storage
    Real *sa = new Real[n];

    // work space
    Real *sf = new Real[n];
    memcpy(sf, sy, n*sizeof(Real));

    // apply Newton's finite difference method
    sa[0] = sf[0];
    for (mkIndex j=1; j<n; j++) {
        for (mkIndex k=n-1; k>=j; k--) {
            Real d = sx[k]-sx[k-j];
            if (d == 0.0)
                sf[k] = 0.0;  // assume that the derivative is 0.0 and apply the Hermite interpolation
            else
                sf[k] = (sf[k]-sf[k-1]) / d;
        }
        sa[j] = sf[j];
    }

    // free work space
    delete [] sf;

    // return the result
    return Vector(n, sa);
}

// return the evaluated P(z0), i.e. the value of P(z) at z=z0, where P(z) is a Newton polynomial defined by
//    P(z) = a(1) + a(2)*(z-x(1)) + a(3)*(z-x(1))*(z-x(2)) + ... + a(n)*(z-x(1))*...*(z-x(n-1))
// this routine also works for evaluating the function value of a Hermite interpolating polynomial,
//    which is in the same form as a Newton polynomial
Real NewtonPolynomialEvaluation(const Vector &a, const Vector &x, const Real z0) {
    // polynomial degree is n-1; number of coefficients is n
    mkIndex n = a.Length();

    // a special case
    if (n == 0u)
        return 0.0;

    // exception handling
    // note that x(1,...,n-1) will do; x(n) is not required
    if (x.Length() < n-1u) {
        *Basics::err << "NewtonPolynomialEvaluation(const Vector&, const Vector&): the number of sampling points is less than the polynomial degree!" << endl;
        Basics::quit(1);
    }

    // evaluate P(z0)
    const Real *sa = a.Store()+n-1;
    const Real *sx = x.Store()+n-2;
    mkIndex jj = n-1;
    Real fval = *sa--;
    while (jj--) {
        fval *= (z0-(*sx--));
        fval += (*sa--);
    }

    // return the result
    return fval;
}

// compute a base filter P(z) which is a continuous, piecewise polynomial P(z) expanded
// in a basis of `translated' (i.e. scale-and-shift) Chebyshev polynomials in each interval
//
// The base filter P(z) equals P_j(z) for z in the j-th interval [intv(j), intv(j+1)), where
// P_j(z) a Hermite interpolating polynomial
//
// input:
// intv is a vector which defines the intervals; the j-th interval is [intv(j), intv(j+1))
// HiLowFlags determines the shape of the base filter P(z)
// Consider the j-th interval [intv(j), intv(j+1)]
// HighLowFlag[j-1]==1,  P(z)==1 for z in [intv(j), intv(j+1)]
//                 ==0,  P(z)==0 for z in [intv(j), intv(j+1)]
//                 ==-1, [intv(j), intv(j+1)] is a transition interval;
//                       P(intv(j)) and P(intv(j+1)) are defined such that P(z) is continuous
// baseDeg is the degree of smoothness of the Hermite (piecewise polynomial) interpolation
// to be precise, the i-th derivative of P(z) is zero, i.e. d^{i}P(z)/dz^i==0, at all interval end points z=intv(j) for i=1,...,baseDeg
//
// output:
// P(z) expanded in a basis of `translated' (scale-and-shift) Chebyshev polynomials
// to be precise, for z in the j-th interval [intv(j),intv(j+1)), P(z) equals
//     P_j(z) = pp(1,j)*S_0(z) + pp(2,j)*S_1(z) + ... + pp(n,j)*S_{n-1}(z),
// where S_i(z) is the `translated' Chebyshev polynomial in that interval,
//     S_i((z-c)/h) = T_i(z),  c = (intv(j)+intv(j+1))) / 2,  h = (intv(j+1)-intv(j)) / 2,
// with T_i(z) the Chebyshev polynomial of the first kind,
//     T_0(z) = 1, T_1(z) = z, and T_i(z) = 2*z*T_{i-1}(z) - T_{i-2}(z) for i>=2
// the return matrix is the matrix of Chebyshev coefficients pp just described
//
// note that the degree of P(z) in each interval is (at most) 2*baseDeg+1, with 2*baseDeg+2 coefficients
// let n be the length of intv; then there are n-1 intervals
// therefore the return matrix pp is of size (2*baseDeg+2)-by-(n-1)
Matrix HermiteBaseFilterInChebyshevBasis(const Vector &intv, const int *HighLowFlags, mkIndex baseDeg) {
    // number of points
    mkIndex npoints = intv.Length();

    if (npoints <= 1) {
        // an exception: no intervals, return an empty matrix
        return Matrix();
    }

    // jj is initialized as the number of intervals
    mkIndex jj = npoints-1u;

    // baseFilter is the pointer to the storage of the output matrix to be computed
    Real *baseFilter = new Real[(2u*baseDeg+2u)*jj];

    // a pointer traverses and fills basefilter[] with the Chebyshev coefficients
    Real *sbf = baseFilter;

    // the main loop to compute the Chebyshev coefficients
    const int *hilo = HighLowFlags;
    const Real *currentPoint = intv.Store();
    while (jj--) {
        // get the flag of the current interval
        Real flag = (Real)(*hilo++);
        if (flag == -1.0) {
            // flag == -1 means that the current interval is a transition polynomial

            // get flag2, the flag of the next interval
            Real flag2 = (Real)(*hilo);
            // the flag of the previous interval is 1-flag2
            Real flag0 = 1.0-flag2;

            // allocate memory for two vectors
            Real *x = new Real[2u*baseDeg+2u];
            Real *y = new Real[2u*baseDeg+2u];

            // two pointers for traversing x[] and y[]
            Real *sx = x;
            Real *sy = y;

            // find the current interval [aa,bb]
            Real aa = *currentPoint++;
            Real bb = *currentPoint;

            // now left-hand side
            mkIndex ii = baseDeg+1u;
            while (ii--) {
                *sy++ = flag0;
                *sx++ = aa;
            }

            // now right-hand side
            ii = baseDeg+1;
            while (ii--) {
                *sy++ = flag2;
                *sx++ = bb;
            }

            // build a Newton polynomial (indeed, the generalized Hermite interpolating polynomial) with x[] and y[]
            Vector px(2u*baseDeg+2u, x);
            Vector pp = NewtonPolynomial(px, Vector(2u*baseDeg+2u, y));

            // pp contains coefficients of the Newton polynomial P(z) in the current interval [aa,bb], where
            // P(z) = pp(1) + pp(2)*(z-px(1)) + pp(3)*(z-px(1))*(z-px(2)) + ... + pp(n)*(z-px(1))*...*(z-px(n-1))

            // translate the Newton coefficients to the Chebyshev coefficients
            Vector qq = ExpandNewtonPolynomialInChebyshevBasis(aa, bb, pp, px);
            // qq contains coefficients of the polynomial in [aa,bb] in the `translated' Chebyshev basis

            // copy the Chebyshev coefficients to baseFilter
            // OCTAVE/MATLAB: B(:,j) = qq, where j = (npoints-1)-jj and B is the return matrix
            const Real *sq = qq.Store();
            ii = 2u*baseDeg + 2u;
            while (ii--)
                *sbf++ = *sq++;
        }
        else {
            // a constant polynomial P(z)=flag, where either flag==0 or flag==1
            // OCTAVE/MATLAB: B(1,j) = flag, where j = (npoints-1)-jj and B is the return matrix
            *sbf ++ = flag;

            // the other coefficients are all zero, since P(z) is a constant
            // OCTAVE/MATLAB: B(1,j) = 0, where j = (npoints-1)-jj and B is the return matrix
            mkIndex ii = 2u*baseDeg+1u;
            while (ii--)
                *sbf++ = 0.0;

            // for the next point
            currentPoint++;
        }
    }

    // return the result
    return Matrix(2u*baseDeg+2u, npoints-1u, baseFilter);
}


////////////////////////////////////////////////////////////////////////////
//    Base Filter
////////////////////////////////////////////////////////////////////////////

// same as the next GetIntervals(), with default IntervalOptions
PolynomialFilterInfo GetIntervals(Vector &intervals, const Vector &frame, mkIndex polyDeg, mkIndex baseDeg) {
    IntervalOptions opts;
    return GetIntervals(intervals, frame, polyDeg, baseDeg, opts);
}

// this routine determines the intervals (including the transition one(s)) by an iterative process
//
// frame is a vector consisting of 4 ordered elements:
//     [frame(1),frame(4)] is the interval which (tightly) contains all eigenvalues, and
//     [frame(2),frame(3)] is the interval in which the eigenvalues are sought
// baseDeg is the left-and-right degree of the base filter for each interval
// polyDeg is the (maximum possible) degree of s(z), with z*s(z) the polynomial filter
// intv is the output; the j-th interval is [intv(j),intv(j+1))
// opts is a collection of interval options
//
// the base filter P(z) is a piecewise polynomial from Hermite interpolation with degree baseDeg at each end point of intervals
//
// the polynomial filter Q(z) is in the form z*s(z), i.e. Q(0)==0, such that ||(1-P(z))-z*s(z)||_w is minimized with
// s(z) a polynomial of degree up to polyDeg
//
// the resulting polynomial filter Q(z) satisfies Q(x)>=Q(y) for x in [frame(2),frame(3)], and
// y in [frame(1),frame(4)] but not in [frame(2),frame(3)]
//
// the routine returns an instance of PolynomialFilterInfo which gives some information of the polynomial filter
PolynomialFilterInfo GetIntervals(Vector &intervals, const Vector &frame, mkIndex polyDeg, mkIndex baseDeg, IntervalOptions &opts) {
    if (frame.Length() != 4) {
        *Basics::err << "GetIntervals(Vector&, const Vector&, mkIndex, mkIndex, IntervalOptions&): the frame vector, the second argument, should contain exactly 4 elements!" << endl;
        Basics::quit(1);
    }
    const Real a = frame(1), a1 = frame(2), b1 = frame(3), b = frame(4);
    // [a,b]   is the range of all eigenvalues
    // [a1,b1] is the range of wanted eigenvalues
    if (a>a1 || a1>b1 || b1>b) {
        *Basics::err << "GetIntervals(Vector&, const Vector&, mkIndex, mkIndex, IntervalOptions&): values in the frame vector, the second argument, should be non-decreasing!" << endl;
        Basics::quit(1);
    }
    if (a1 == b1) {
        *Basics::err << "GetIntervals(Vector&, const Vector&, mkIndex, mkIndex, IntervalOptions&): the range of wanted eigenvalues cannot be of size zero!" << endl;
        Basics::quit(1);
    }
    PolynomialFilterInfo filterInfo;
    filterInfo.filterType = 2;      // mid-pass filter, for interior eigenvalues
    if (b == b1) {
        if (a == a1) {
            *Basics::err << "GetIntervals(Vector&, const Vector&, mkIndex, mkIndex, IntervalOptions&): a polynomial filter should not cover all eigenvalues!" << endl;
            Basics::quit(1);
        }
        filterInfo.filterType = 1;  // high-pass filter, for largest eigenvalues
    }
    else if (a == a1) {
        // filterType = 3, low-pass filter, for smallest eigenvalues
        // however it should have been converted to the largest eigenvalue problem
        *Basics::err << "GetIntervals(Vector&, const Vector&, mkIndex, mkIndex, IntervalOptions&): filterType==3 for smallest eigenvalues should be pre-converted to filterType==1 for largest eigenvalues!" << endl;
        Basics::quit(1);
    }

    // the following recipe follows Yousef Saad (2005, 2006) with a few minor adaptations / enhancements
    // initialization
    static const int HighLowFlags[] = { 1, -1, 0, -1, 1 };  // if filterType is 1, only first 3 elements 1,-1,0 will be used
    Real halfPlateau = 0.5*(b1-a1)*opts.initialPlateau;     // half length of the "plateau" of the (dual) base filter
    Real leftDelta = (b1-a1)*opts.initialShiftStep;         // initial left shift
    Real rightDelta = leftDelta;                            // initial right shift
    opts.numGridPoints = Basics::max(opts.numGridPoints, (mkIndex)(2.0*(b-a)/halfPlateau));
    Real gridSize = (b-a) / (Real)(opts.numGridPoints);
    Vector intv;
    if (filterInfo.filterType == 2) {
        // for interior eigenvalues
        intv.resize(6);
        intv(1) = a;
        intv(6) = b;
        // intv(2), intv(3), intv(4), and intv(5) to be determined
    }
    else {
        // filterType == 1 (or 3 with conversion), for extreme eigenvalues
        intv.resize(4);
        intv(1) = a;
        intv(4) = b;
        // intv(2) and intv(3) to be determined
    }
    Matrix baseFilter, polyFilter;
    mkIndex numIter, numLeftSteps = 1, numRightSteps = 1;
    const mkIndex numLookMore = 2u*(mkIndex)(0.5+(log(2.0)/log(opts.shiftStepExpansionRate)));
    mkIndex numMoreLooked = 0;
    Real yLimit, ySummit;
    Real yLimitGap = 0.0, yLeftSummit = 0.0, yLeftBottom = 0.0, yRightSummit = 0.0, yRightBottom = 0.0;  // initialization is dispensable
    Real z1 = a1 - leftDelta;
    Real z2 = b1 + rightDelta;
    filterInfo.filterOK = 0;  // not yet found any OK filter
    bool leftBoundaryMet = false, rightBoundaryMet = false;

    // initialize the intervals, mainly for the case opts.maxOuterIter == 0
    intervals = intv;
    intervals(2) = z1;
    if (filterInfo.filterType == 2) {
        // a mid-pass filter for interior eigenvalues
        intervals(5) = z2;
        Real c = (a1+b1) / 2.0;
        intervals(3) = c - halfPlateau;
        intervals(4) = c + halfPlateau;
    }
    else {
        // filterType == 1 (or 3 with conversion) for extreme eigenvalues
        intervals(3) = z1 + (b1-z1)*opts.transitionIntervalRatio;
    }

    // the main loop
    for (numIter=1; numIter<=opts.maxOuterIter; numIter++) {
        // outer loop updates z1 and z2
        if (z1 <= a) {
            z1 = a;
            leftBoundaryMet = true;
        }
        // a <= z1 < (a1)
        if (filterInfo.filterType == 2) {
            // a mid-pass filter for interior eigenvalues
            if (z2 >= b) {
                z2 = b;
                rightBoundaryMet = true;
            }
            Real c = (z1+z2) / 2.0;
            // a <= z1 < c-h < c+h < z2 <= b, where h is halfPlateau
            // [z1, c-h] and [c+h, z2] are transition interval
            intv(2) = z1;
            intv(5) = z2;
            Real c1 = z1 + halfPlateau;
            intv(3) = z1;           // i.e. c1 - halfPlateau
            intv(4) = c1 + halfPlateau;
            baseFilter = HermiteBaseFilterInChebyshevBasis(intv, HighLowFlags, baseDeg);
            polyFilter = FilteredConjugateResidualPolynomial(baseFilter, intv, opts.intervalWeights, polyDeg);
            // Real fc1 = PiecewisePolynomialEvaluationInChebyshevBasis(polyFilter, intv, b1) - PiecewisePolynomialEvaluationInChebyshevBasis(polyFilter, intv, a1);
            Real c2 = z2 - halfPlateau;
            intv(3) = c2 - halfPlateau;
            intv(4) = z2;           // i.e. c2 + halfPlateau
            baseFilter = HermiteBaseFilterInChebyshevBasis(intv, HighLowFlags, baseDeg);
            polyFilter = FilteredConjugateResidualPolynomial(baseFilter, intv, opts.intervalWeights, polyDeg);
            Real fc2 = PiecewisePolynomialEvaluationInChebyshevBasis(polyFilter, intv, b1) - PiecewisePolynomialEvaluationInChebyshevBasis(polyFilter, intv, a1);
            yLimitGap = Basics::infinity;
            mkIndex ii = opts.maxInnerIter;
            while (ii-- && !(yLimitGap <= opts.yLimitTol)) {
                // recursive bisection to get c such that polynomial values at a1 and b1 are balanced; i.e. P(a1) are P(b1) approximately the same
                c = (c1+c2) / 2.0;
                intv(3) = c - halfPlateau;
                intv(4) = c + halfPlateau;
                baseFilter = HermiteBaseFilterInChebyshevBasis(intv, HighLowFlags, baseDeg);
                polyFilter = FilteredConjugateResidualPolynomial(baseFilter, intv, opts.intervalWeights, polyDeg);
                Real fc = PiecewisePolynomialEvaluationInChebyshevBasis(polyFilter, intv, b1) - PiecewisePolynomialEvaluationInChebyshevBasis(polyFilter, intv, a1);
                if (fc*fc2 < 0.0) {
                    c1 = c;
                    // fc1 = fc;
                }
                else {
                    c2 = c;
                    fc2 = fc;
                }
                yLimitGap = Basics::abs(fc);
            }
        }
        else {  // filterType == 1 (or 3 with conversion) for extreme eigenvalues
            intv(2) = z1;
            intv(3) = z1 + (b1-z1)*opts.transitionIntervalRatio;
            baseFilter = HermiteBaseFilterInChebyshevBasis(intv, HighLowFlags, baseDeg);
            polyFilter = FilteredConjugateResidualPolynomial(baseFilter, intv, opts.intervalWeights, polyDeg);
        }
        // polyFilter contains the coefficients of the polynomial filter which approximates phi(x) expanded in the `translated' Chebyshev basis
        // psi(x) = 1.0 - phi(x) is the dual base filter approximated by a polynomial in the form x*p(x)
        Real yLeftLimit  = 1.0 - PiecewisePolynomialEvaluationInChebyshevBasis(polyFilter, intv, a1);
        Real yRightLimit = 1.0 - PiecewisePolynomialEvaluationInChebyshevBasis(polyFilter, intv, b1);
        yLimit  = (yLeftLimit < yRightLimit) ? yLeftLimit : yRightLimit;
        ySummit = (yLeftLimit > yRightLimit) ? yLeftLimit : yRightLimit;
        Real x = a1, y;
        while ((x+=gridSize) < b1) {
            y = 1.0 - PiecewisePolynomialEvaluationInChebyshevBasis(polyFilter, intv, x);
            if (y < yLimit)
                yLimit = y;
            if (y > ySummit)
                ySummit = y;
        }
        // now yLimit is the minimum of x*p(x) for x in [a1, b1]
        bool stepLeft = false;
        bool stepRight = false;
        if ((yLimit < yLeftLimit && yLimit < yRightLimit) || yLimit < opts.yBottomLine) {
            // very bad, step to see what will happen
            stepLeft = true;
            if (filterInfo.filterType == 2)
                stepRight = true;
        }
        else if (filterInfo.filterType == 2) {
            if (yLeftLimit < yRightLimit) {
                if (yRightLimit-yLeftLimit > opts.yLimitTol)
                    stepLeft = true;
            }
            else if (yLeftLimit-yRightLimit > opts.yLimitTol)
                stepRight = true;
        }
        if (!stepLeft) {
            yLeftBottom = yLeftLimit;
            x = a1;
            while ((x-=gridSize) >= a) {
                y = 1.0 - PiecewisePolynomialEvaluationInChebyshevBasis(polyFilter, intv, x);
                if (y < yLeftBottom)
                    yLeftBottom = y;
                else if (y > yLeftBottom)
                    break;
            }
            yLeftSummit = yLeftBottom;
            while ((x-=gridSize) >= a) {
                y = 1.0 - PiecewisePolynomialEvaluationInChebyshevBasis(polyFilter, intv, x);
                if (y > yLeftSummit) {
                    yLeftSummit = y;
                    if (yLeftSummit > yLimit*opts.yRippleLimit) {
                        stepLeft = true;
                        break;
                    }
                }
                if (y < yLeftBottom)
                    yLeftBottom = y;
            }
        }
        if (filterInfo.filterType == 2 && !stepRight) {
            yRightBottom = yRightLimit;
            x = b1;
            while ((x+=gridSize) <= b) {
                y = 1.0 - PiecewisePolynomialEvaluationInChebyshevBasis(polyFilter, intv, x);
                if (y < yRightBottom)
                    yRightBottom = y;
                else if (y > yRightBottom)
                    break;
            }
            yRightSummit = yRightBottom;
            while ((x+=gridSize) <= b) {
                y = 1.0 - PiecewisePolynomialEvaluationInChebyshevBasis(polyFilter, intv, x);
                if (y > yRightSummit) {
                    yRightSummit = y;
                    if (yRightSummit > yLimit*opts.yRippleLimit) {
                        stepRight = true;
                        break;
                    }
                }
                if (y < yRightBottom)
                    yRightBottom = y;
            }
        }
        if (!stepLeft && !stepRight) {
            Real bottom;
            if (filterInfo.filterType == 2)  // a mid-pass filter
                bottom = Basics::min(yLeftBottom, yRightBottom);
            else  // a high-pass filter (or a low-pass filter w/ conversion); there is no right wing
                bottom = yLeftBottom;
            Real qIndex = 1.0 - (yLimit-bottom) / (ySummit-bottom);
            if (filterInfo.filterOK == 0 || filterInfo.filterQualityIndex < qIndex) {
                // found the first OK filter or a better filter
                intervals = intv;
                filterInfo.filterOK = 1;
                filterInfo.filterQualityIndex = qIndex;
                filterInfo.numIter = numIter;
                filterInfo.yLimit = yLimit;
                filterInfo.ySummit = ySummit;
                filterInfo.numLeftSteps = numLeftSteps;
                filterInfo.yLeftSummit = yLeftSummit;
                filterInfo.yLeftBottom = yLeftBottom;
                if (filterInfo.filterType == 2) {
                    filterInfo.yLimitGap = yLimitGap;
                    filterInfo.numRightSteps= numRightSteps;
                    filterInfo.yRightSummit = yRightSummit;
                    filterInfo.yRightBottom = yRightBottom;
                }
                numMoreLooked = 0;
            }
            else if (++numMoreLooked == numLookMore) {
                // filter has been optimal
                filterInfo.filterOK = 2;
                break;
            }
            // try stepping further to see whether it can improve
            stepLeft = true;
            if (filterInfo.filterType == 2)
                stepRight = true;
        }
        // check whether we can really "step"
        if (leftBoundaryMet) {
            if (filterInfo.filterType == 1 || rightBoundaryMet)
                break;  // cannot step further, so break the loop
            if (stepLeft) {
                // cannot step left, so try stepping right
                stepLeft = false;
                stepRight = true;
            }
        }
        if (rightBoundaryMet && stepRight) {
            // cannot step right, so try stepping left
            stepRight = false;
            stepLeft = true;
        }
        // now "step"
        if (stepLeft) {
            numLeftSteps ++;
            if (filterInfo.filterType == 2)
                leftDelta *= opts.shiftStepExpansionRate;  // expand the moving step for faster convergence
            z1 -= leftDelta;
        }
        if (stepRight) {
            numRightSteps ++;
            rightDelta *= opts.shiftStepExpansionRate;     // expand the shift step for faster convergence
            z2 += rightDelta;
        }
        if (filterInfo.filterType == 2) {
            // shrink the "plateau" of the (dual) base filter
            if (stepLeft && stepRight)
                halfPlateau /= opts.plateauShrinkRate;
            else
                halfPlateau /= sqrt(opts.plateauShrinkRate);
        }
    }
    filterInfo.totalNumIter = numIter;
    return filterInfo;
}


////////////////////////////////////////////////////////////////////////////
//    Chebyshev Polynomials
////////////////////////////////////////////////////////////////////////////

// translate the coefficients of a Newton polynomial to the coefficients
// in a basis of the `translated' (scale-and-shift) Chebyshev polynomials
//
// input:
// a Newton polynomial defined by vectors a and x:
//     P(z) = a(1) + a(2)*(z-x(1)) + a(3)*(z-x(1))*(z-x(2)) + ... + a(n)*(z-x(1))*...*(z-x(n-1))
// the interval [aa,bb] defines the `translated' Chebyshev polynomials S_i(z) = T_i((z-c)/h),
//     where c=(aa+bb)/2 and h=(bb-aa)/2, and T_i is the Chebyshev polynomial of the first kind
// note that T_i is defined by T_0(z)=1, T_1(z)=z, and T_i(z)=2*z*T_{i-1}(z)+T_{i-2}(z) for i>=2
//
// output:
// a vector q containing the Chebyshev coefficients:
//     P(z) = q(1)*S_0(z) + q(2)*S_1(z) + ... + q(n)*S_{n-1}(z)
Vector ExpandNewtonPolynomialInChebyshevBasis(Real aa, Real bb, const Vector &a, const Vector &x) {
    // number of coefficients
    mkIndex n = a.Length();

    // exception handling
    if (x.Length()+1u < n) {
        // x(n) is not required, but we need x(1,...,n-1)
        *Basics::err << "ExpandNewtonPolynomialInChebyshevBasis(Real, Real, const Vector&, const Vector&): the number of sampling points is less than the polynomial degree!" << endl;
        Basics::quit(1);
    }
    if (n == 0)
        return Vector();

    // pointers for traversing a and x
    const Real *sa = a.Store()+n;
    const Real *sx = x.Store()+n-1;

    // storage for output
    Real *q = new Real[n];

    // work space, will be freed
    Real *q2 = new Real[n];

    // set q[0] = a(n)
    *q = *--sa;

    // the main loop for translation
    Real c = (aa+bb)/2.0;
    Real h = (bb-aa)/2.0;
    Real h2 = h/2.0;
    for (mkIndex m=1; m<=n-1; m++) {
        // compute q2[0:m-1] = (c-x[n-m-1])*q[0:m-1]
        mkIndex mm = m;
        Real *sq = q;
        Real *sq2 = q2;
        Real c2 = c-(*--sx);
        while (mm--)
            *(sq2++) = c2 * (*sq++);
        *sq2 = 0.0;         // q2[m] = 0.0
        *(q2+1) += h*(*q);  // q2[1] = q2[1] + h*q[0]

        // compute q2[0:m-2] = q2[0:m-2] + h2*q[1:m-1]
        mm = m-1;
        sq2 = q2;
        sq = q+1;
        while (mm--)
            *(sq2++) += h2 * (*sq++);

        // compute q2[2:m] = q2[2:m] + h2*q[1:m-1]
        mm = m-1;
        sq2 = q2+2;
        sq = q+1;
        while (mm--)
            *(sq2++) += h2*(*sq++);

        // compute q[0:m] = q2[0:m]
        mm = m+1;
        sq2 = q2;
        sq = q;
        while (mm--)
            *sq++ = *sq2++;
        *q += (*--sa);      // q[0] = q[0] + p[n-m-1]
    }

    // free work space
    delete [] q2;

    // return the result
    return Vector(n, q);
}

// evaluate P(z) at z=z0, where P(z) is a polynomial expanded in a basis of
// the `translated' (i.e. scale-and-shift) Chebyshev polynomials
//
// input:
// c is a vector of Chebyshev coefficients which defines the polynomial
//     P(z) = c(1)*S_0(z) + c(2)*S_1(z) + ... + c(n)*S_{n-1}(z),
// where S_i is the `translated' Chebyshev polynomial S_i((z-c)/h) = T_i(z), with
//     c = (intv(j)+intv(j+1)) / 2,  h = (intv(j+1)-intv(j)) / 2
// note that T_i(z) is the Chebyshev polynomial of the first kind,
//     T_0(z) = 1, T_1(z) = z, and T_i(z) = 2*z*T_{i-1}(z) - T_{i-2}(z) for i>=2
//
// output:
// the evaluated value of P(z) at z=z0
Real PolynomialEvaluationInChebyshevBasis(const Vector &c, Real z0, Real aa, Real bb) {
    // number of coefficients, one more than polynomial degree
    mkIndex deg1 = c.Length();

    // a special: zero polynomial
    if (deg1 == 0u)
        return 0.0;

    // find the translated z0
    Real zz;
    if (aa==-1.0 && bb==1.0)
        zz = z0;  // treat it as a special case to reduce rounding errors
    else
        zz = (aa==bb) ? 0.0 : -1.0+2.0*(z0-aa)/(bb-aa);

    // compute y = P(z0), where we utilize the Chebyshev recursion
    const Real *sc = c.Store();
    Real y = *sc++;
    Real t0, t1 = 1.0, t2 = zz;
    mkIndex ii = deg1-1u;
    while (ii--) {
        // Chebyshev recursion: T_0(zz)=1, T_1(zz)=zz, and T_{i+1}(zz) = 2*zz*T_i(zz) + T_{i-1}(zz) for i>=2
        // the values of T_{i+1}(zz), T_i(zz), T_{i-1}(zz) are stored in t0, t1, t2, respectively
        t0 = 2*zz*t1 - t2;
        // it also works for the base case / the first iteration, where t0 equals 2*zz*1-zz == zz which is T_1(zz)
        t2 = t1;
        t1 = t0;
        y += (*sc++) * t0;
    }

    // return the evaluated P(z0)
    return y;
}

// evaluate P(z) at z=z0, where P(z) is a piecewise polynomial expanded
// in a basis of the (optionally translated, i.e. scale-and-shift) Chebyshev polynomials for each interval
//
// input:
// intv is a vector which defines the intervals; the j-th interval is [intv(j), intv(j+1))
// pp is a matrix of Chebyshev coefficients which defines a piecewise polynomial P(z)
// in a basis of the `translated' Chebyshev polynomials in each interval
// the polynomial P_j(z) in the j-th interval, i.e. when z is in [intv(j), intv(j+1)), is defined by the j-th column of pp:
//     if basisTranslated == false, then
//         P_j(z) = pp(1,j)*T_0(z) + pp(2,j)*T_1(z) + ... + pp(n,j)*T_{n-1}(z),
//     where T_i(z) is the Chebyshev polynomial of the first kind,
//         T_0(z) = 1, T_1(z) = z, and T_i(z) = 2*z*T_{i-1}(z) - T_{i-2}(z) for i>=2
//     if basisTranslated == true, then
//         P_j(z) = pp(1,j)*S_0(z) + pp(2,j)*S_1(z) + ... + pp(n,j)*S_{n-1}(z),
//     where S_i is the `translated' Chebyshev polynomial S_i((z-c)/h) = T_i(z), with
//         c = (intv(j)+intv(j+1)) / 2,  h = (intv(j+1)-intv(j)) / 2
//
// output:
// the evaluated value of P(z) at z=z0
//
// note that if z0 falls below the first interval, then the polynomial in the first interval will be used for evaluating P(z0)
//           if z0 flies over  the last  interval, then the polynomial in the last  interval will be used for evaluating P(z0)
Real PiecewisePolynomialEvaluationInChebyshevBasis(const Matrix &pp, const Vector &intv, Real z0, bool basisTranslated) {
    // number of points, one more than number of intervals
    mkIndex npoints = intv.Length();

    // determine the interval which contains z0
    const Real *sintv = intv.Store()+1;
    mkIndex idx = 1;
    if (npoints > 2 && z0 >= *sintv) {
        sintv ++;
        while (++idx < npoints-1) {
            if (z0 < *sintv)
                break;
            sintv ++;
        }
    }
    // idx==1 if npoints<=2; otherwise idx satisfies:
    //     intv(idx) <= z0 < intv(idx+1),  if 2 <= idx <= npoints-2
    //     z0 < intv(idx+1),               if idx == 1
    //     intv(idx) <= z0,                if idx == npoints-1
    // in addition, sintv points to &intv(idx+1)

    // exception handling
    if (pp.Ncols() < idx) {
        *Basics::err << "PiecewisePolynomialEvaluationInChebyshevBasis(const Matrix&, const Vector&, Real, bool): fewer polynomials (columns of the first argument) than the intervals (length of the second argument)!" << endl;
        Basics::quit(1);
    }

    if (basisTranslated) {
        // the basis consists of `translated' Chebyshev polynomials
        // find the interval of concern, [aa,bb]
        Real aa = *(sintv-1), bb = *sintv;
        return PolynomialEvaluationInChebyshevBasis(pp.column(idx), z0, aa, bb);
    }
    // else the basis consists of standard Chebyshev polynomials, with interval [-1.0,1.0] for integration
        return PolynomialEvaluationInChebyshevBasis(pp.column(idx), z0);
}

// compute the weighted inner product of two piecewise polynomials expanded
// in a basis of `translated' (i.e. scale-and-shift) Chebyshev polynomials for each interval
//
// pp and qq are two matrices of Chebyshev coefficients which define the piecewise polynomials P(z) and Q(z), respectively
// for z in the j-th interval, P(z) equals
//     P_j(z) = pp(1,j)*S_0(z) + pp(2,j)*S_1(z) + ... + pp(n,j)*S_{n-1}(z),
// and Q(z) equals
//     Q_j(z) = qq(1,j)*S_0(z) + qq(2,j)*S_1(z) + ... + qq(n,j)*S_{n-1}(z),
// where S_i(z) is the `translated' Chebyshev polynomial in that interval,
//     S_i((z-c)/h) = T_i(z),  c = (aa+bb)) / 2,  h = (bb-aa) / 2,
// with T_i(z) the Chebyshev polynomial of the first kind
//     T_0(z) = 1, T_1(z) = z, and T_i(z) = 2*z*T_{i-1}(z) - T_{i-2}(z) for i>=2
//
// the (scaled) j-th interval inner product is defined by
//     <P_j,Q_j> = (Pi/2)*( pp(1,j)*qq(1,j) + sum_{k} pp(k,j)*qq(k,j) ),
// which comes from the property
//     <T_0,T_0>=pi, <T_i,T_i>=pi/2 for i>=1, and <T_i,T_j>=0 for i!=j
//
// the weighted inner product is <P,Q> = sum_{j} intervalWeights(j)*<P_j,Q_j>,
// which is the return value
//
// note that for unit weights, pass an empty vector of intervalWeights (i.e. of length 0)
Real PiecewisePolynomialInnerProductInChebyshevBasis(const Matrix &pp, const Matrix &qq, const Vector &intervalWeights) {
    // number of effective coefficients, one more than the effective polynomial degree
    mkIndex deg1 = Basics::min(pp.Nrows(), qq.Nrows());

    // a special case, zero polynomial
    if (deg1 == 0)
        return 0.0;

    // number of intervals
    mkIndex numIntv = Basics::min(pp.Ncols(), qq.Ncols());

    // the extra amount to skip, if any
    mkIndex dp = pp.Nrows() - deg1;
    mkIndex dq = qq.Nrows() - deg1;

    // pointers for traversing p and q
    const Real *sp = pp.Store();
    const Real *sq = qq.Store();

    Real ans = 0.0;
    if (intervalWeights.Length() == 0) {
        // unit weights
        // compute ans = sum_{i=1,...,numIntv} [ pp(1,i)*qq(1,i) + sum_{k=1,...,deg} pp(k,i)*qq(k,i) ]
        mkIndex ii = numIntv;
        while (ii--) {
            // sum up  pp(1,i)*qq(1,i) + sum_{k=1,...,deg} pp(k,i)*qq(k,i), where i = numIntv-ii
            ans += (*sp) * (*sq);  // the first term pp(1,i)*qq(1,i) is being added twice
            mkIndex kk = deg1;
            while (kk--) {
                // add pp(k,i)*qq(k,i), where k = deg1-kk
                ans += (*sp++) * (*sq++);
            }
            // skip the extra
            sp += dp;
            sq += dq;
        }
    }
    else if (intervalWeights.Length() < numIntv) {
        *Basics::err << "PiecewisePolynomialInnerProductInChebyshevBasis(const Matrix&, const Matrix&, const Vector&): there are fewer weights than intervals!" << endl;
        Basics::quit(1);
    }
    else {
        // scaled by intervalWeights(i) in the i-th interval
        // compute ans = sum_{i=1,...,numIntv} intervalWeights(i)*[ pp(1,i)*qq(1,i) + sum_{k=1,...,deg} pp(k,i)*qq(k,i) ]
        mkIndex ii = numIntv;
        const Real *sw = intervalWeights.Store();
        while (ii--) {
            // compute ans2 = pp(1,i)*qq(1,i) + sum_{k=1,...,deg} pp(k,i)*qq(k,i), where i = numIntv-ii
            Real ans2 = (*sp) * (*sq);  // the first term pp(1,i)*qq(1,i) is being added twice
            mkIndex kk = deg1;
            while (kk--) {
                // add pp(k,i)*qq(k,i), where k = deg1-kk
                ans2 += (*sp++) * (*sq++);
            }
            // compute ans += ans2*intervalWeights(i)
            ans += ans2*(*sw++);
            // skip the extra
            sp += dp;
            sq += dq;
        }
    }

    // return the inner product
    return ans*Basics::pi/2.0;
}

// compute Q(z) = z*P(z), where P(z) and Q(z) are piecewise polynomials expanded
// in a basis of `translated' (i.e. scale-and-shift) Chebyshev polynomials for each interval
//
// P(z) and Q(z) are stored as matrices of Chebyshev coefficients pp and qq, respectively
//
// For z in the j-th interval, P(z) equals
//     P_j(z) = pp(1,j)*S_0(z) + pp(2,j)*S_1(z) + ... + pp(n,j)*S_{n-1}(z),
// and Q(z) equals
//     Q_j(z) = qq(1,j)*S_0(z) + qq(2,j)*S_1(z) + ... + qq(n,j)*S_{n-1}(z),
// where S_i(z) is the `translated' Chebyshev polynomial in that interval,
//     S_i((z-c)/h) = T_i(z),  c = (intv(j)+intv(j+1))) / 2,  h = (intv(j+1)-intv(j)) / 2,
// with T_i(z) the Chebyshev polynomial of the first kind
//     T_0(z) = 1, T_1(z) = z, and T_i(z) = 2*z*T_{i-1}(z) - T_{i-2}(z) for i>=2
//
// the returned matrix is qq which represents Q(z) = z*P(z)
Matrix PiecewisePolynomialInChebyshevBasisMultiplyX(const Matrix &pp, const Vector &intv) {
    // We exploit the relation x*T_i(x) = (T_{i+1}(x) + T_{i-1}(x))/2 and S_i(z) = T_i((z-c)/h).
    // Let x=(z-c)/h. We obtain z*S_i(z) = 0.5*h*(S_{i-1}(z)+S_{i+1}(z)) + c*S_i(z), for i>=1.
    // For i=0, x*T_i(x) = T_{i+1}(x). Therefore, z*S_i(z) = h*S_{i+1}(z) + c*S_i(z).

    // number of intervals
    mkIndex nintv = pp.Ncols();

    // exception handling
    if (intv.Length() != nintv+1u) {
        *Basics::err << "PiecewisePolynomialInChebyshevBasisMultiplyX(const Matrix &pp, const Vector &intv): number of (piecewise) polynomials pp.Ncols() must agree with number of intervals intv.Length()-1!" << endl;
        Basics::quit(1);
    }

    // number of coefficients, one more than the polynomial degree
    mkIndex deg1 = pp.Nrows();

    // a special case, zero polynomial
    if (deg1 == 0) {
        // zero polynomial times z is still zero polynomial
        return Matrix();
    }

    // straightforward implementation, not best efficient
    /*
    Matrix qq(deg1+1, nintv);
    for (mkIndex jj=1; jj<=nintv; jj++) {
        Real c = (intv(jj) + intv(jj+1)) / 2.0;
        Real h = (intv(jj+1) - intv(jj)) / 2.0;
        Real h2 = h / 2.0;
        for (mkIndex ii=1; ii<=deg1; ii++)
            qq(ii,jj) = c*pp(ii,jj);
        qq(deg1+1,jj) = 0.0;
        mkIndex ii = 1;
        qq(ii+1,jj) += h*pp(ii,jj);
        for (mkIndex ii=2; ii<=deg1; ii++) {
            qq(ii-1,jj) += h2*pp(ii,jj);
            qq(ii+1,jj) += h2*pp(ii,jj);
        }
    }
    return qq;
    */
    // the following implementation achieves better efficiency

    // pointers for traversing pp and intervals
    const Real *sp = pp.Store();
    const Real *sp2 = pp.Store();
    const Real *sintv = intv.Store();

    // allocate memory
    Real *q = new Real[(deg1+1u)*nintv];
    Real *sq = q;
    Real *sq2 = q+1;

    mkIndex jj = nintv;
    while (jj--) {
        // consider interval between intv(j) and intv(j+1), where j == nintv-jj

        // compute c = (intv(j) + intv(j+1))/2
        Real c = 0.5*(*sintv + *(sintv+1));

        // compute h = (intv(j+1) - intv(j))/2  and  h2 = h/2
        Real h = 0.5*(*(sintv+1) - (*sintv));
        Real h2 = 0.5*h;

        // compute q(1:deg1,j) = c*p(1:deg1,j)
        mkIndex i = deg1;
        while (i--)
            *sq++ = c*(*sp++);

        // set q(deg1+1,j) = 0.0
        *sq++ = 0.0;

        // compute q(2,j) = q(2,j) + h*p(1,j)
        *(sq2++) += h*(*sp2++);

        // compute q(3:deg1+1,j) = q(3:deg1+1,j) + h2*p(2:deg1,j) and
        //    then q(1:deg1-1,j) = q(1:deg1-1,j) + h2*p(2:deg1,j)
        i = deg1-1;
        while (i--) {
            Real tmp = h2*(*sp2++);
            *(sq2-2) += tmp;
            *(sq2++) += tmp;
        }

        // for pointing to q(2,j+1)
        sq2 ++;
        // for the next interval
        sintv ++;
    }

    // return the result
    return Matrix(deg1+1u, nintv, q);
}


////////////////////////////////////////////////////////////////////////////
//    Conjugate Residual Method in the Polynomial Space
////////////////////////////////////////////////////////////////////////////

// this routine employs a conjugate-residual-type algorithm in polynomial space to minimize ||P(z)-Q(z)||_w,
// where P(z), the base filter, is the input piecewise polynomial, and
//       Q(z) is the output polynomial satisfying Q(0)==1, i.e. the constant term of Q(z) is 1
// niter is the number of conjugate-residual iterations; therefore, the degree of Q(z) is up to niter+1
// both P(z) and Q(z) are expanded in the `translated' (scale-and-shift) Chebyshev basis for each interval,
// and presented as matrices of Chebyshev coefficients, denoted by pp and qq, respectively
//
// input:
// intv is a vector which defines the intervals; the j-th interval is [intv(j),intv(j+1))
// w is a vector of Chebyshev weights; the weight of j-th interval is w(j)
//     the interval weights define the inner product of two continuous functions and then
//     the derived w-norm ||P(z)-Q(z)||_w
// pp is a matrix of Chebyshev coefficients which defines the piecewise polynomial P(z)
// to be specific, for z in [intv(j), intv(j+1)), P(z) equals
//     P_j(z) = pp(1,j)*S_0(z) + pp(2,j)*S_1(z) + ... + pp(niter+2,j)*S_{niter+1}(z),
// where S_i(z) is the `translated' Chebyshev polynomial in that interval,
//     S_i((z-c)/h) = T_i(z),  c = (intv(j)+intv(j+1))) / 2,  h = (intv(j+1)-intv(j)) / 2,
// with T_i(z) the Chebyshev polynomial of the first kind,
//     T_0(z) = 1, T_1(z) = z, and T_i(z) = 2*z*T_{i-1}(z) - T_{i-2}(z) for i>=2
//
// output:
// the return matrix, denoted by qq, represents a polynomial Q(z) with degree up to 1+niter
// and satisfying Q(0)==1, such that ||P(z))-Q(z)||_w is minimized
// this polynomial Q(z) is expanded in the `translated' Chebyshev basis for each interval
// to be precise, considering z in [intv(j), intv(j+1)), Q(z) equals
//     Q_j(z) = qq(1,j)*S_0(z) + qq(2,j)*S_1(z) + ... + qq(niter+2,j)*S_{niter+1}(z)
//
// note:
// 1. since Q(0)==1, P(0)==1 is expected; if P(0)!=1, one can translate P(z)
//    for example, if P(0)==0, one can use 1-P(z) as input instead of P(z)
// 2. typically, the base filter, defined by pp and intv, is from Hermite interpolation
//    in intervals [intv(j),intv(j+1)) for j=1,...,nintv, with nintv the number of intervals
Matrix FilteredConjugateResidualPolynomial(const Matrix &baseFilter, const Vector &intv, const Vector &intervalWeights, mkIndex niter) {
    if (intv.Length() <= 1u) {
        // exception, no intervals, no polynomial
        return Matrix();
    }

    // initialize polynomial ppol to be 1 (i.e. multiplicative identity) in all intervals
    // number of intervals
    mkIndex nintv = intv.Length() - 1u;
    // allocate memory
    Real *pp = new Real[2u*nintv];
    // set the polynomial be 1 for all intervals
    Real *sp = pp;
    mkIndex jj = nintv;
    while (jj--) {
        *sp++ = 1.0;
        *sp++ = 0.0;
    }
    Matrix ppol(2, nintv, pp);  // the initial p-polynomial (corresponding to the A-conjugate vector p in CG)

    Matrix rpol = ppol;         // rpol is the r-polynomial (corresponding to the residual vector r in CG)
    Matrix cpol = ppol;         // cpol is the "corrected" residual polynomial
    Matrix appol = PiecewisePolynomialInChebyshevBasisMultiplyX(ppol, intv);
    Matrix arpol = appol;
    Real rho = PiecewisePolynomialInnerProductInChebyshevBasis(rpol, arpol, intervalWeights);
    for (mkIndex i=0; i<niter; i++) {
        Real den = PiecewisePolynomialInnerProductInChebyshevBasis(appol, appol, intervalWeights);
        Real alp0 = rho / den;
        Real alp = alp0 - PiecewisePolynomialInnerProductInChebyshevBasis(baseFilter, appol, intervalWeights) / den;
        rpol = Matrix::xsum(rpol, (-alp0)*appol);
        cpol = Matrix::xsum(cpol, (-alp)*appol);
        if (i+1u == niter)
            break;
        arpol = PiecewisePolynomialInChebyshevBasisMultiplyX(rpol, intv);
        Real rho0 = rho;
        rho = PiecewisePolynomialInnerProductInChebyshevBasis(rpol, arpol, intervalWeights);
        Real bet = rho / rho0;
        ppol = Matrix::xsum(rpol, bet*ppol);
        appol = Matrix::xsum(arpol, bet*appol);
    }

    // return the (corrected) residual polynomial
    return cpol;
}

// this routine employs a conjugate-residual-type algorithm in polynomial space to compute
// x = x0 + s(A)*r0 with r0 = b - A*x0, such that ||(1-P(z))-z*s(z)||_w is minimized, where
// P(z) is a given piecewise polynomial, called the base filter,
// s(z) is a polynomial of degree up to niter, the number of conjugate-residual iterations,
// and b and x0 are given vectors
//
// note that P(z) is expanded in the `translated' (scale-and-shift) Chebyshev basis for each interval,
// and presented as a matrix of Chebyshev coefficients pp
//
// input:
// A is a sparse matrix
// x0, b are vectors
// niter is the number of conjugate-residual iterations
// intv is a vector which defines the intervals; the j-th interval is [intv(j),intv(j+1))
// w is a vector of Chebyshev weights; the weight of j-th interval is w(j)
//     the interval weights define the inner product of two continuous functions and then
//     the derived w-norm ||P(z)-Q(z)||_w
// pp is a matrix of Chebyshev coefficients which defines the piecewise polynomial P(z)
// to be specific, for z in [intv(j), intv(j+1)), P(z) equals
//     P_j(z) = pp(1,j)*S_0(z) + pp(2,j)*S_1(z) + ... + pp(niter+2,j)*S_{niter+1}(z),
// where S_i(z) is the `translated' Chebyshev polynomial in that interval,
//     S_i((z-c)/h) = T_i(z),  c = (intv(j)+intv(j+1))) / 2,  h = (intv(j+1)-intv(j)) / 2,
// with T_i(z) the Chebyshev polynomial of the first kind,
//     T_0(z) = 1, T_1(z) = z, and T_i(z) = 2*z*T_{i-1}(z) - T_{i-2}(z) for i>=2
// tol is the tolerance; if the residual polynomial in z-norm is dropped by a factor lower
//     than tol, then stop the conjugate-residual iteration
//
// output:
// the return vector is x = x0 + s(A)*r0 with r0 = b - A*x0, such that ||(1-P(z))-z*s(z)||_w is minimized,
// subject to that s(z) is a polynomial of degree up to niter, where P(z) is the base filter
// in short, z*s(z) approximates 1-P(z)
//
// note:
// 1. since z*s(z) approximates 1-P(z), P(0)==1 is expected; if P(0)!=1, one can translate P(z)
//    for example, if P(0)==0, one can use 1-P(z) as input instead of P(z)
// 2. typically, the base filter, defined by pp and intv, is from Hermite interpolation
//    in intervals [intv(j),intv(j+1)) for j=1,...,nintv, with nintv the number of intervals
// 3. a common application is to compute R(A)*b, where R(z) approximates 1-P(z)
//    in this case, one can set x0 = 0 and then the return vector is x = s(A)*b, where
//    z*s(z) approximates 1-P(z); therefore, A*x is the wanted R(A)*b
Vector FilteredConjugateResidualMatrixPolynomialVectorProduct(const SparseMatrix &A, const Vector &x0, const Vector &b,
                                                 const Matrix &baseFilter, const Vector &intv, const Vector &intervalWeights,
                                                 mkIndex niter, Real tol) {
    if (intv.Length() <= 1u) {
        // exception, no intervals, no update on x0
        return x0;
    }

    // initialize polynomial ppol to be 1 (i.e. multiplicative identity) in all intervals
    // number of intervals
    mkIndex nintv = intv.Length()-1u;
    // allocate memory
    Real *pp = new Real[2u*nintv];
    // set the polynomial be 1 for all intervals
    Real *sp = pp;
    mkIndex jj = nintv;
    while (jj--) {
        *sp++ = 1.0;
        *sp++ = 0.0;
    }
    Matrix ppol(2, nintv, pp);  // the initial p-polynomial (corresponding to the A-conjugate vector p in CG)

    // corrected CR in polynomial space
    Matrix rpol = ppol;         // rpol is the r-polynomial (corresponding to the residual vector r in CG)
    Matrix cpol = ppol;         // cpol is the "corrected" residual polynomial
    Matrix appol = PiecewisePolynomialInChebyshevBasisMultiplyX(ppol, intv);
    Matrix arpol = appol;
    Real rho00 = PiecewisePolynomialInnerProductInChebyshevBasis(rpol, arpol, intervalWeights);  // initial rho
    Real rho = rho00;

    // corrected CR in vector space
    Vector x = x0;
    #ifdef USE_CSRMV
        Vector r = b-(x*A);
    #else
        Vector r = b-(A*x);
    #endif
    Vector p = r;
    #ifdef USE_CSRMV
        Vector ap = p*A;
    #else
        Vector ap = A*p;
    #endif
    // alp0, alp, and bet will be from the polynomial space in the iteration

    // to trace the residual errors, execute the lines starting with "//*"
    //* Vector err(iter+1);
    //* err(1) = rho;
    for (mkIndex i=0; i<niter; i++) {
        // iteration in the polynomial space
        Real den = PiecewisePolynomialInnerProductInChebyshevBasis(appol, appol, intervalWeights);
        Real alp0 = rho / den;
        Real alp = alp0 - PiecewisePolynomialInnerProductInChebyshevBasis(baseFilter, appol, intervalWeights) / den;
        rpol = Matrix::xsum(rpol, (-alp0)*appol);
        cpol = Matrix::xsum(cpol, (-alp)*appol);
        arpol = PiecewisePolynomialInChebyshevBasisMultiplyX(rpol, intv);
        Real rho0 = rho;
        rho = PiecewisePolynomialInnerProductInChebyshevBasis(rpol, arpol, intervalWeights);

        // update x in the vector space
        x += alp*p;
        //* err(i+2) = (b-A*x).Norm2();
        if (rho < tol*rho00)
            break;

        // finish the iteration in the polynomial space
        Real bet = rho / rho0;
        ppol = Matrix::xsum(rpol, bet*ppol);
        appol = Matrix::xsum(arpol, bet*appol);

        // finish the iteration in the vector space
        r -= alp0*ap;
        p = r + bet*p;
        #ifdef USE_CSRMV
            ap = r*A + bet*ap;  // the only matrix-vector product in the loop
        #else
            ap = A*r + bet*ap;  // the only matrix-vector product in the loop
        #endif
    }
    return x;
}

// variables for polynomial filter interface
// a pointer to the input sparse symmetric matrix via setFilter()
const SparseMatrix *PolynomialFilterInterface::Sptr = NULL;
// the sparse symmetric matrix S in computing R(S)*v; S is translated from the input sparse symmetric matrix S0 in setFilter()
SparseMatrix PolynomialFilterInterface::S;
// the j-th interval is [intervals(j),intervals(j+1))
// the intervals are decided by GetIntervals() invoked by setFilter() which takes parameters frame, baseDeg, and polyDeg:
// 1. frame is a vector consisting of 4 ordered elements:
//     [frame(1),frame(4)] is the interval which (tightly) contains all eigenvalues of S, and
//     [frame(2),frame(3)] is the interval in which the eigenvalues are sought
// 2. baseDeg is the left-and-right degree of the base filter for each interval
// 3. polyDeg is the (maximum possible) degree of z*s(z), with s(z) the polynomial filter
Vector PolynomialFilterInterface::intervals;
// intervalWeights(j) is the weight of j-th interval [intervals(j), intervals(j+1))
Vector PolynomialFilterInterface::intervalWeights;
// maximum possible degree of s(z), with z*s(z) the polynomial filter
mkIndex PolynomialFilterInterface::polyDegree;
// left-and-right degree of the base filter in each interval
mkIndex PolynomialFilterInterface::baseDegree;
// a base filter, which is a piecewise polynomial typically from Hermite interpolation
// for j-th interval, the polynomial is expanded in the `translated' (scale-and-shift) Chebyshev basis,
// with the Chebyshev coefficients stored in baseFilter(:,j)
Matrix PolynomialFilterInterface::baseFilter;
// a zero vector of the length the number of rows/columns of S
Vector PolynomialFilterInterface::zn;
// the information of computing the polynomial filter; see class PolynomialFilterInfo in polyfilt.h for more information
PolynomialFilterInfo PolynomialFilterInterface::filterInfo;

// this member function sets the sparse matrix S translated from S0, and also the base filter P(z)
// S0 is the sparse symmetric matrix of concern
// frame is a vector of 4 ordered elements:
//     [frame(1),frame(4)] is the interval which (tightly) contains all eigenvalues of S0, and
//     [frame(2),frame(3)] is the interval in which the eigenvalues of S0 are sought
// polyDeg is the degree of s(z), with z*s(z) the polynomial filter
// baseDeg is the left and right degrees of the base filter for each interval
// opts is a collection of options to determine the intervals
Vector PolynomialFilterInterface::setFilter(const SparseMatrix &S0, const Vector &frame, mkIndex polyDeg, mkIndex baseDeg, IntervalOptions &opts) {
    if (frame.Length() != 4) {
        *Basics::err << "PolynomialFilterInterface::setFilter(const Matrix&, const Vector&, mkIndex, mkIndex): the frame vector, the second argument, must contain 4 values (for 3 intervals)!" << endl;
        Basics::quit(1);
    }
    intervalWeights = opts.intervalWeights;
    baseDegree = baseDeg;
    polyDegree = polyDeg;
    mkIndex n = S0.Nrows();
    Vector intervals2;
    if (frame(1) == frame(2)) {
        // low pass filter, convert it to high pass filter
        S = frame(4)*SparseMatrix::eye(n) - S0;
        Sptr = &S;
        Vector frame2 = frame(4) - frame.reverse();
        filterInfo = GetIntervals(intervals, frame2, polyDegree, baseDegree, opts);
        // translate the interval back for return
        intervals2 = frame(4) - intervals.reverse();
    }
    else {
        // it can be a mid-pass filter or a high-pass filter
        if (frame(1) == 0.0) {
            Sptr = &S0;
            filterInfo = GetIntervals(intervals, frame, polyDegree, baseDegree, opts);
            // not translation of intervals
            intervals2 = intervals;
        }
        else {
            S = S0 - frame(1)*SparseMatrix::eye(n);
            Sptr = &S;
            Vector frame2 = frame - frame(1);
            filterInfo = GetIntervals(intervals, frame2, polyDegree, baseDegree, opts);
            // shift the intervals back for return
            intervals2 = intervals + frame(1);
        }
    }
    static const int HighLowFlags[] = { 1, -1, 0, -1, 1 };
    baseFilter = HermiteBaseFilterInChebyshevBasis(intervals, HighLowFlags, baseDegree);

    // set the vector of zeros
    zn = Vector::zeros(n);

    return intervals2;
}
