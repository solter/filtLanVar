#ifndef POLYFILT_H
#define POLYFILT_H
/*! \file polyfilt.h
 *  Routines for polynomial filtering used in projection-type methods for the symmetric / Hermitian eigenvalue problems.
 */

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


////////////////////////////////////////////////////////////////////////////////
//    Newton - Hermite Polynomial Interpolation
////////////////////////////////////////////////////////////////////////////////
/** @defgroup group_poly_intp Polynomial Interpolation
 *
 *  @{
 */

//! Build a polynomial <em>P</em>(<em>z</em>) by Newton's divided differences and return the coefficient vector <em>a</em>.
/*! A Newton polynomial is in the form
 *  <em>P</em>(<em>z</em>) = <em>a</em>(<em>1</em>) + <em>a</em>(<em>2</em>)*(<em>z-x</em>(<em>1</em>))
 *                         + <em>a</em>(<em>3</em>)*(<em>z-x</em>(<em>1</em>))*(<em>z-x</em>(<em>2</em>))<em> + ... + a</em>(<em>n</em>)*(<em>z-x</em>(<em>1</em>))*...*(<em>z-x</em>(<em>n-1</em>)),
 *  such that <em>P</em>(<em>x</em>(<em>i</em>))<em> = y</em>(<em>i</em>) for <em>i=1,...,n</em>.
 *  \param  x         is the input vector of the domain values.
 *  \param  y         is the input vector of the range values.
 *
 *  \return  The vector <em>a</em> of coefficients.
 *
 *  All vectors <em>a</em>, <em>x</em>, and <em>y</em> are of the same length <em>n</em>.
 *  \remark  If <em>x</em>(<em>i</em>)==<em>x</em>(<em>j</em>) for some <em>i\f$\neq\f$j</em>, then
 *           it is assumed that the derivative of <em>P</em>(<em>z</em>) is to be zero at <em>x</em>(<em>i</em>)
 *           and the Hermite polynomial interpolation is applied.
 *  \remark  In general, if there are <em>k</em> <em>x</em>(<em>i</em>)'s with the same value <em>x0</em>, then
 *           the <em>j</em>th order derivative of <em>P</em>(<em>z</em>) is zero at <em>z=x0</em> for <em>j=1,...,k-1</em>.
 */
// build P(z) by Newton's divided differences in the form
//     P(z) = a(1) + a(2)*(z-x(1)) + a(3)*(z-x(1))*(z-x(2)) + ... + a(n)*(z-x(1))*...*(z-x(n-1)),
// such that P(x(i)) = y(i) for i=1,...,n, where
//     x,y are input vectors of length n, and a is the output vector of length n
// if x(i)==x(j) for some i!=j, then it is assumed that the derivative of P(z) is to be zero at x(i),
//     and the Hermite polynomial interpolation is applied
// in general, if there are k x(i)'s with the same value x0, then
//     the j-th order derivative of P(z) is zero at z=x0 for j=1,...,k-1
Vector NewtonPolynomial(const Vector &x, const Vector &y);

//! Evaluate <em>P</em>(<em>z0</em>), i.e. the value of <em>P</em>(<em>z</em>) at <em>z=z0</em>, where <em>P</em>(<em>z</em>) is a Newton polynomial defined by <em>a</em> and <em>x</em>.
/*! The newton polynomial is
 *  <em>P</em>(<em>z</em>) = <em>a</em>(<em>1</em>) + <em>a</em>(<em>2</em>)*(<em>z-x</em>(<em>1</em>))
 *                         + <em>a</em>(<em>3</em>)*(<em>z-x</em>(<em>1</em>))*(<em>z-x</em>(<em>2</em>))<em> + ... + a</em>(<em>n</em>)*(<em>z-x</em>(<em>1</em>))*...*(<em>z-x</em>(<em>n-1</em>)).
 *
 *  \remark  This routine also works for evluating the function value of a Hermite interpolating polynomial,
 *           which is in the same form as a Newton polynomial.
 *
 *  \return  The evaluated <em>P</em>(<em>z0</em>), i.e. the value of <em>P</em>(<em>z</em>) at <em>z=z0</em>.
 */
// return the evaluated P(z0), i.e. the value of P(z) at z=z0, where P(z) is a Newton polynomial defined by
//    P(z) = a(1) + a(2)*(z-x(1)) + a(3)*(z-x(1))*(z-x(2)) + ... + a(n)*(z-x(1))*...*(z-x(n-1))
// this routine also works for evluating the function value of a Hermite interpolating polynomial,
//    which is in the same form as a Newton polynomial
Real NewtonPolynomialEvaluation(const Vector &a, const Vector &x, const Real z0);

/** @} */  // end of group_poly_intp



////////////////////////////////////////////////////////////////////////////
//    Base Filter
////////////////////////////////////////////////////////////////////////////
/** @defgroup group_base_filter Base Filter
 *
 *  @{
 */
// begin of group_base_filter
//! Compute a base filter <em>P</em>(<em>z</em>) which is a continuous, piecewise polynomial <em>P</em>(<em>z</em>) expanded in a basis of `translated' (i.e. scale-and-shift) Chebyshev polynomials in each interval.
/*! The base filter <em>P</em>(<em>z</em>) equals <em>P<sub>j</sub></em>(<em>z</em>) for <em>z</em> in the <em>j</em>th interval [<em>intv</em>(<em>j</em>), <em>intv</em>(<em>j+1</em>)), where
 *  <em>P<sub>j</sub></em>(<em>z</em>) a Hermite interpolating polynomial.
 *
 *  \param  intv        is a vector which defines the intervals. The <em>j</em>th interval is [<em>intv</em>(<em>j</em>), <em>intv</em>(<em>j+1</em>)).
 *  \param  HiLowFlags  determines the shape of the base filter P(z). Consider the <em>j</em>th interval.
 *          <ul>
 *          <li>  If <em>HighLowFlag</em>[<em>j-1</em>]==<em>1</em>,  then <em>P</em>(<em>z</em>)==<em>1</em> for <em>z</em> in [<em>intv</em>(<em>j</em>), <em>intv</em>(<em>j+1</em>)].
 *          <li>  If <em>HighLowFlag</em>[<em>j-1</em>]==<em>1</em>,  then <em>P</em>(<em>z</em>)==<em>1</em> for <em>z</em> in [<em>intv</em>(<em>j</em>), <em>intv</em>(<em>j+1</em>)].
 *          <li>  If <em>HighLowFlag</em>[<em>j-1</em>]==<em>-1</em>, then [<em>intv</em>(<em>j</em>), <em>intv</em>(<em>j+1</em>)] is a transition interval;
 *                                        <em>q</em>(<em>intv</em>(<em>j</em>)) and <em>q</em>(<em>intv</em>(<em>j+1</em>)) are defined such that <em>q</em>(<em>z</em>) is continuous.
 *          </ul>
 *  \param  baseDeg      defines the smoothness of the piecewise polynomial <em>P</em>(<em>z</em>). To be precise,
 *          To be precise, the <em>i</em>th derivative of <em>P</em>(<em>z</em>) is zero, i.e. <em>d<sup>i</sup>P(z)/dz<sub>i</sub>==0</em>,
 *          at all interval end points <em>z=intv</em>(<em>j</em>) for <em>i=1,...,baseDeg</em>.
 *
 *  The base filter <em>P</em>(<em>z</em>) expanded in a basis of `translated' (scale-and-shift) Chebyshev polynomials,
 *  presented by a matrix of Chebyshev coefficients <em>pp</em>.
 *
 *  For <em>z</em> in the <em>j</em>th interval [<em>intv</em>(<em>j</em>), <em>intv</em>(<em>j+1</em>)), <em>P</em>(<em>z</em>) equals
 *  <em>P<sub>j</sub></em>(<em>z</em>) = <em>pp</em>(<em>1,j</em>)*<em>S<sub>0</sub></em>(<em>z</em>) + <em>pp</em>(<em>2,j</em>)*<em>S<sub>1</sub></em>(<em>z</em>)<em> + ... + pp</em>(<em>n,j</em>)<em>*S<sub>n-1</sub></em>(<em>z</em>).
 *
 *  Here <em>S<sub>i</sub></em>(<em>z</em>) is the `translated' Chebyshev polynomial in the <em>j</em>th interval, with
 *         <em>S<sub>i</sub></em>((<em>z-c</em>)/<em>h</em>) = <em>T<sub>i</sub></em>(<em>z</em>),  <em>c</em> = (<em>intv</em>(<em>j</em>) + <em>intv</em>(<em>j+1</em>))/<em>2</em>,  <em>h</em> = (<em><em>intv</em>(<em>j+1</em>) - <em>intv</em>(<em>j</em>)</em>)/<em>2</em>,\n
 *  where <em>T<sub>i</sub></em>(<em>z</em>) is the Chebyshev polynomial of the first kind,
 *         <em>T<sub>0</sub></em>(<em>z</em>) = <em>1</em>, <em>T<sub>1</sub></em>(<em>z</em>) = <em>z</em>, and <em>T<sub>i</sub></em>(<em>z</em>) = <em>2*z*T<sub>i-1</sub></em>(<em>z</em>) - <em>T<sub>i-2</sub></em>(<em>Z</em>) for <em>i>=2</em>.
 *
 *  \return The matrix <em>pp</em> of Chebyshev coefficients.
 *
 *  \remark The degree of <em>P</em>(<em>z</em>) in each interval is (at most) <em>2*baseDeg+1</em>, with <em>2*baseDeg+2</em> coefficients.
 *  \remark Let <em>n</em> be the length of <em>intv</em>; then there are <em>n-1</em> intervals.
 *          The return matrix <em>pp</em> is of size (<em>2*baseDeg+2</em>)-by-(<em>n-1</em>).
 */
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
// to be preicse, for z in the j-th interval [intv(j),intv(j+1)), P(z) equals
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
Matrix HermiteBaseFilterInChebyshevBasis(const Vector &intv, const int *HiLowFlags, mkIndex baseDeg);

//! An instance of this class, taken by GetIntervals() and PolynomialFilterInterface::setFilter(), is a collection of options to determine the intervals which decides the base filter.
// an instance of this class, taken by GetIntervals() and PolynomialFilterInterface::setFilter(), is a collection of options to determine the intervals which decides the base filter
// it is used as an input argument of GetIntervals(), and also
struct IntervalOptions {
    //! Interval weights (default <em>100,1,1,1,100</em>).
    /*! \remark  For mid-pass filters, there are <em>5</em> intervals, including <em>2</em> transition intervals.
     *           So all the <em>5</em> weights <em>intervalWeights</em>(<em>1</em>),...,<em>intervalWeights</em>(<em>5</em>) are used.
     *  \remark  For high-pass filters (and low-pass filters with conversion), there are <em>3</em> intervals, including <em>1</em> transition interval.
     *           So only the first <em>3</em> weights <em>intervalWeights</em>(<em>1</em>),<em>intervalWeights</em>(<em>2</em>),<em>intervalWeights</em>(<em>3</em>) are used.
     */
    // interval weights, default 100,1,1,1,100
    Vector intervalWeights;

    // the following parameter is for high-pass filters (and low-pass filters with conversion) only
    //! The (relative) length of transition interval (default <em>0.6</em>), effective only for high-pass filters, and low-pass filters with conversion.
    // the (relative) length of transition interval (default 0.6)
    Real transitionIntervalRatio;

    // the following parameters are for mid-pass filters only
    //! This parameter tells whether to reverse the interval or not. It is effective only for mid-pass filters (default <b>false</b>).
    // this parameter tells whether to reverse the interval or not; it is effective only for mid-pass filters (default false)
    bool reverseInterval;
    //! Initial length of `plateau', relative to the length of the interval of desired eigenvalues (default <em>0.1</em>), effective only for mid-pass filters.
    // initial length of `plateau' relative to the length of the interval of desired eigenvalues (default 0.1)
    Real initialPlateau;
    //! The rate at which the plateau shrinks at each iteration (default <em>1.5</em>), effective only for mid-pass filters.
    // the rate at which the plateau shrinks at each iteration (default 1.5)
    Real plateauShrinkRate;
    //! Initial shift step, relative to the length of the interval of desired eigenvalues (default <em>0.01</em>), effective only for mid-pass filters.
    // initial shift step, relative to the length of the interval of desired eigenvalues (default 0.01)
    Real initialShiftStep;
    //! The rate at which the shift step expands (default <em>1.5</em>), effective only for mid-pass filters.
    // the rate at which the shift step expands (default 1.5)
    Real shiftStepExpansionRate;
    //! Maximum number of inner iterations (default <em>30</em>), effective only for mid-pass filters.
    // maximum number of inner iterations (default 30)
    mkIndex maxInnerIter;
    //! A mid-pass filter <em>P</em>(<em>x</em>) should have <em>P</em>(<em>a1</em>)<em>=P</em>(<em>b1</em>), where [<em>a1,b1</em>] is the interval of desired eigenvalues;\n <em>yLimitTol</em> (default <em>0.0001</em>) is the tolerance allowed for |<em>P</em>(<em>a1</em>)-<em>P</em>(<em>b1</em>)|.
    // a mid-pass filter P(x) should have P(a1)=P(b1), where [a1,b1] is the interval of desired eigenvalues
    // yLimitTol (default 0.0001) is the tolerance allowed for |P(a1)-P(b1)|
    Real yLimitTol;

    // the following parameters are for all filters
    //! Maximum number of outer iterations (default <em>50</em>).
    // maximum number of outer iterations, default 50
    mkIndex maxOuterIter;
    //! Number of grid points, used to measure the maximum function value of the polynomial <em>P</em>(<em>z</em>) for <em>z</em> not in the interval in which the eigenvalues are sought (default <em>200</em>).
    /*! \remark  The value of <em>numGridPoints</em> will automatically be increased if it is too small.
     */
    // number of grid points, used to measure the maximum function value of the polynomial P(z) for z not in the interval which the eigenvalues are sought (default 200)
    // the value of numGridPoints will automatically be increased if it is too small
    mkIndex numGridPoints;
    //! This is the bottom line (default <em>0.001</em>) which the function value of the polynomial must be greater than for <em>z</em> is in the interval of desired eigenvalues.
    // this is the bottom line (default 0.001) which the function value of the polynomial must be greater than for z is in the interval of desired eigenvalues
    Real yBottomLine;
    //! The limit of height of ripples (not in the interval of desired eigenvalues) relative to the bottom of polynomial values for <em>z</em> in the interval of desired eigenvalues (default <em>0.8</em>).
    // the limit of height of ripples (not in the interval of desired eigenvalues) relative to the bottom of polynomial values for z in the interval of desired eigenvalues (default 0.8)
    Real yRippleLimit;

    //! A default constructor to set default parameters.
    // a default constructor to set default parameters
    IntervalOptions() {
        // set default interval weights
        intervalWeights.resize(5);
        intervalWeights(1) = 100.0;
        intervalWeights(2) = 1.0;
        intervalWeights(3) = 1.0;
        intervalWeights(4) = 1.0;
        intervalWeights(5) = 100.0;

        // for high-pass filters (and low-pass filters with conversion)
        transitionIntervalRatio = 0.6;

        // for mid-pass filters
        reverseInterval = false;
        initialPlateau = 0.1;
        plateauShrinkRate = 1.5;
        initialShiftStep = 0.01;
        shiftStepExpansionRate = 1.5;
        maxOuterIter = 50;
        yLimitTol = 0.0001;

        // for all filters
        maxInnerIter = 30;
        numGridPoints = 200;
        yBottomLine = 0.001;
        yRippleLimit = 0.8;
    }
};

//! The routine GetIntervals() returns an instance of this class, which gives the information of a polynomial filter.
// the routine GetIntervals() returns an instance of this class, which gives the information of a polynomial filter
struct PolynomialFilterInfo {
    // the following returned values are for all filters
    //! The partition of the range of spectrum which decides the base filter and the polynomial filter.
    // The partition of the range of spectrum which decides the base filter and the polynomial filter.
    Vector intervals;
    //! <em>1</em> means a high-pass filter (or low-pass filter with conversion);\n <em>2</em> means a mid-pass filter.
    // 1 means a high-pass filter (or low-pass filter with conversion); 2 means a mid-pass filter
    int filterType;
    //! <em>0</em> means no acceptable is found; <em>1</em> means an OK filter is found; <em>2</em> means an optimal filter is found.
    // 0 means no acceptable is found; 1 means an OK filter is found; 2 means an optimal filter is found
    int filterOK;
    //! Between <em>0.0</em> and <em>1.0</em>; the higher the better quality of the filter.
    // between 0.0 and 1.0; the higher the better quality of the filter
    Real filterQualityIndex;
    //! Number of iterations to get the (transition) intervals.
    // number of iterations to get the (transition) intervals
    mkIndex numIter;
    //! Total number of iterations performed.
    /*! This may include a few extra iterations without improving the (transition) intervals, compared to PolynomialFilterInfo::numIter.
     */
    // total number of iterations performed; this may include a few extra iterations without improving the (transition) intervals, compared to numIter
    mkIndex totalNumIter;
    //! The lowest polynomial value <em>P</em>(<em>z</em>) for <em>z</em> in the interval [<em>a1</em>,<em>b1</em>] of desired eigenvalues.
    // the lowest polynomial value P(z) for z in the interval [a1,b1] of desired eigenvalues
    Real yLimit;
    //! The height of (highest, if more than one) summit in the interval [<em>a1</em>,<em>b1</em>] of desired eigenvalues.
    // the height of (highest, if more than one) summit in the interval [a1,b1] of desired eigenvalues
    Real ySummit;
    //! Number of steps moving leftward.
    // number of steps moving leftward
    mkIndex numLeftSteps;
    //! The height of highest summit in the left-hand side of the interval of desired eigenvalues.
    // the height of highest summit in the left-hand side of the interval of desired eigenvalues
    Real yLeftSummit;
    //! The height of lowest bottom in the left-hand side of the interval of desired eigenvalues.
    // the height of lowest bottom in the left-hand side of the interval of desired eigenvalues
    Real yLeftBottom;

    // the following returned values are for high-pass filters (or low-pass filters after conversion)
    //! In general, it is |<em>p</em>(<em>a1</em>)-<em>p</em>(<em>b1</em>)|, where [<em>a1</em>,<em>b1</em>] is the interval of desired eigenvalues.
    // in general it is |p(a1)-p(b1)|, where [a1,b1] is the interval of desired eigenvalues
    Real yLimitGap;
    //! Number of steps moving rightward.
    // number of steps moving rightward
    mkIndex numRightSteps;
    //! The height of highest summit in the right-hand side of the interval of desired eigenvalues.
    // the height of highest summit in the right-hand side of the interval of desired eigenvalues
    Real yRightSummit;
    //! The height of lowest bottom in the right-hand side of the interval of desired eigenvalues.
    // the height of lowest bottom in the right-hand side of the interval of desired eigenvalues
    Real yRightBottom;

    //! A default constructor to initialze all zero.
    // a default constructor to initialize all zero
    PolynomialFilterInfo() {
        // zero initialization
        filterOK = 0;
        filterQualityIndex = 0.0;

        numIter = 0;
        totalNumIter = 0;
        yLimit = 0.0;
        yLimitGap = 0.0;
        ySummit = 0.0;

        numLeftSteps = 0;
        yLeftSummit = 0.0;
        yLeftBottom = 0.0;

        numRightSteps = 0;
        yRightSummit = 0.0;
        yRightBottom = 0.0;
    }
};

//! This routine determines the intervals (including the transition one(s)) by an interative process.
/*! \param  frame     is a vector consisting of 4 ordered elements:\n
 *          <ul>
 *          <li>  [<em>frame</em>(<em>1</em>),<em>frame</em>(<em>4</em>)] is the interval which (tightly) contains all eigenvalues.
 *          <li>  [<em>frame</em>(<em>2</em>),<em>frame</em>(<em>3</em>)] is the interval in which the eigenvalues are sought.
 *          </ul>
 *  \param  baseDeg   is the left-and-right degree of the base filter for each interval.
 *  \param  polyDeg   is the (maximum possible) degree of <em>s</em>(<em>z</em>), with <em>z*s</em>(<em>z</em>) the polynomial filter.
 *  \param  intv      is the output; the <em>j</em>th interval is [<em>intv</em>(<em>j</em>),<em>intv</em>(<em>j+1</em>)).
 *  \param  opts      is a collection of interval options.
 *
 *  The base filter <em>P</em>(<em>z</em>) is a piecewise polynomial from Hermite interpolation with
 *  degree <em>baseDeg</em> at each end point of intervals.
 *
 *  The polynomial filter <em>Q</em>(<em>z</em>) is in the form <em>z*s</em>(<em>z</em>),
 *  i.e. <em>Q</em>(<em>0</em>)<em>==0</em>, such that ||(<em>1-P</em>(<em>z</em>))<em>-z*s</em>(<em>z</em>)||<em><sub>w</sub></em> is minimized with
 *  s(z) a polynomial of degree up to <em>polyDeg</em>.
 *
 *  The resulting polynomial filter <em>Q</em>(<em>z</em>) satisfies
 *  <em>Q</em>(<em>x</em>)>=<em>Q</em>(<em>y</em>)
 *  for <em>x</em> in [<em>frame</em>(<em>2</em>),<em>frame</em>(<em>3</em>)], and
 *  <em>y</em> in [<em>frame</em>(<em>1</em>),<em>frame</em>(<em>4</em>)] but not in [<em>frame</em>(<em>2</em>),<em>frame</em>(<em>3</em>)].
 *
 *  \return An instance of PolynomialFilterInfo which gives some information of the polynomial filter.
 */
// this routine determines the intervals (including the transition one(s)) by an interative process
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
PolynomialFilterInfo GetIntervals(Vector &intv, const Vector &frame, mkIndex polyDeg, mkIndex baseDeg, IntervalOptions &opts);

//! Same as GetIntervals(), with default IntervalOptions.
// same as GetIntervals() above, with default IntervalOptions
PolynomialFilterInfo GetIntervals(Vector &intervals, const Vector &frame, mkIndex polyDeg, mkIndex baseDeg);  // default interval options will be used

/** @} */
// end of group_base_filter



////////////////////////////////////////////////////////////////////////////
//    Chebyshev Polynomials
////////////////////////////////////////////////////////////////////////////

/** @defgroup group_chebyshev Chebyshev Polynomials
 *
 *  @{
 */

//! Translate the coefficients of a Newton polynomial to the coefficients in a basis of the `translated' (scale-and-shift) Chebyshev polynomials.
/*! \param  a,x      define a Newton polynomial as follows: \n
 *          <em>P</em>(<em>z</em>)<em> = a</em>(<em>1</em>)<em> + a</em>(<em>2</em>)*(<em>z-x</em>(<em>1</em>))<em> + a</em>(<em>3</em>)*(<em>z-x</em>(<em>1</em>))*(<em>z-x</em>(<em>2</em>))<em> + ... + a</em>(<em>n</em>)*(<em>z-x</em>(<em>1</em>))*...*(<em>z-x</em>(<em>n-1</em>)).
 *  \param  [aa,bb]  is the interval which defines the `translated' Chebyshev polynomials <em>S<sub>i</sub></em>(<em>z</em>) = <em>T<sub>i</sub></em>((<em>z-c</em>)<em>/h</em>),
 *          where <em>c=(aa+bb)/2</em> and <em>h=(bb-aa)/2</em>, and \n
 *          <em>T<sub>i</sub>(<em>z</em>)</em> is the Chebyshev polynomial of the first kind
 *          <em>T<sub>0</sub></em>(<em>z</em>)=<em>1, T<sub>1</sub></em>(<em>z</em>)<em>=z</em>, and <em>T<sub>i</sub></em>(<em>z</em>)<em>=2*z*T<sub>i-1</sub></em>(<em>z</em>)<em>+T<sub>i-2</sub></em>(<em>Z</em>) for <em>i>=2</em>.
 *
 *  \return A vector <em>q</em> containing the Chebyshev coefficients, such that
 *          <em>P</em>(<em>z</em>) = <em>q</em>(<em>1</em>)<em>*S<sub>0</sub></em>(<em>z</em>)<em> + q</em>(<em>2</em>)<em>*S<sub>1</sub></em>(<em>z</em>)<em> + ... + q</em>(<em>n</em>)<em>*S<sub>n-1</sub></em>(<em>z</em>).
 */
// translate the coefficients of a Newton polynomial to the coefficients
// in a basis of the `translated' Chebyshev polynomials
//
// input:
// a Newton polynomial defined by vectors a and x:
//     P(z) = a(1) + a(2)*(z-x(1)) + a(3)*(z-x(1))*(z-x(2)) + ... + a(n)*(z-x(1))*...*(z-x(n-1))
// the interval [aa,bb] defines the `translated' Chebyshev polynomials S_i(z) = T_i((z-c)/h),
//     where c=(aa+bb)/2 and h=(bb-aa)/2, and T_i is the Chebyshev polynomial of the first kind
// note that T_i is defined by T_0(z)=1, T_1(z)=z, and T_i(z)=2*z*T_{i-1}(z)+T_{i-2}(z) for i>=2
//
// output:
// a vector q containing the Chebyshev coefficients, such that
//     P(z) = q(1)*S_0(z) + q(2)*S_1(z) + ... + q(n)*S_{n-1}(z)
Vector ExpandNewtonPolynomialInChebyshevBasis(Real aa, Real bb, const Vector &a, const Vector &x);

//! Evaluate <em>P</em>(<em>z</em>) at <em>z=z0</em>, where <em>P</em>(<em>z</em>) is a polynomial expanded in a basis of the `translated' (i.e. scale-and-shift) Chebyshev polynomials.
/*! \param  c         is a vector of Chebyshev coefficients which defines the polynomial
 *                    <em>P</em>(<em>z</em>) = <em>c</em>(<em>1</em>)<em>*S<sub>0</sub></em>(<em>z</em>)<em> + c</em>(<em>2</em>)*<em>S<sub>1</sub></em>(<em>z</em>)<em> + ... + c</em>(<em>n</em>)<em>*S<sub>n-1</sub></em>(<em>z</em>),
 *                    where \n
 *                    <em>S<sub>i</sub></em>(<em>z</em>) is the `translated' Chebyshev polynomial <em>S<sub>i</sub></em>((<em>z-c</em>)/<em>h</em>)<em> = T<sub>i</sub></em>(<em>z</em>), with
 *                    <em>c = </em>(<em>aa+bb</em>)<em> / 2</em> and
 *                    <em>h = </em>(<em>bb-aa</em>)<em> / 2</em>.
 *  \param  [aa,bb]   is the interval which defines the `translated' Chebyshev polynomials <em>S<sub>i</sub></em>(<em>z</em>) = <em>T<sub>i</sub></em>((<em>z-c</em>)<em>/h</em>),
 *                    where <em>c=(aa+bb)/2</em> and <em>h=(bb-aa)/2</em>, and \n
 *                    <em>T<sub>i</sub>(<em>z</em>)</em> is the Chebyshev polynomial of the first kind
 *                    <em>T<sub>0</sub></em>(<em>z</em>)=<em>1, T<sub>1</sub></em>(<em>z</em>)<em>=z</em>, and <em>T<sub>i</sub></em>(<em>z</em>)<em>=2*z*T<sub>i-1</sub></em>(<em>z</em>)<em>+T<sub>i-2</sub></em>(<em>Z</em>) for <em>i>=2</em>.
 *  \param  z0        is a number at which the piecewise polynomial <em>P</em>(<em>z</em>) is evaluated.
 *
 *  \return The evaluated value of <em>P</em>(<em>z</em>) at <em>z=z0</em>.
 */
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
Real PolynomialEvaluationInChebyshevBasis(const Vector &c, Real z0, Real aa=-1.0, Real bb=1.0);

//! Evaluate <em>P</em>(<em>z</em>) at <em>z=z0</em>, where <em>P</em>(<em>z</em>) is a piecewise polynomial expanded in a basis of the (optionally translated, i.e. scale-and-shift) Chebyshev polynomials for each interval.
/*! \param  intv      is a vector which defines the intervals. The <em>j</em>th interval is [<em>intv</em>(<em>j</em>), <em>intv</em>(<em>j+1</em>)).
 *  \param  basisTranslated  tells whether the basis of Chebyshev polynomials are translated (<em>basisTranslated==</em><b>true</b>) or not (<em>basisTranslated==</em><b>false</b>).
 *  \param  pp        is a matrix of Chebyshev coefficients which defines a piecewise polynomial in a basis of the (optionally translated) Chebyshev polynomials in each interval. \n
 *          The polynomial <em>P<sub>j</sub></em>(<em>z</em>) in the <em>j</em>th interval, i.e. when <em>z</em> is in [<em>intv</em>(<em>j</em>), <em>intv</em>(<em>j+1</em>)), is defined by the <em>j</em>th column of pp.
 *          <ul>
 *          <li>  If <em>basisTranslated==</em><b>false</b>, then
 *                <em>P<sub>j</sub></em>(<em>z</em>)<em> = pp</em>(<em>1,j</em>)<em>*T<sub>0</sub></em>(<em>z</em>)<em> + pp</em>(<em>2,j</em>)<em>*T<sub>1</sub></em>(<em>z</em>)<em> + ... + pp</em>(<em>n,j</em>)<em>*T<sub>n-1</sub></em>(<em>z</em>), where \n
 *                <em>T<sub>i</sub></em>(<em>z</em>) is the Chebyshev polynomial of the first kind
 *                <em>T<sub>0</sub></em>(<em>z</em>)=<em>1, T<sub>1</sub></em>(<em>z</em>)<em>=z</em>, and <em>T<sub>i</sub></em>(<em>z</em>)<em>=2*z*T<sub>i-1</sub></em>(<em>z</em>)<em>+T<sub>i-2</sub></em>(<em>Z</em>) for <em>i>=2</em>.
 *          <li>  If <em>basisTranslated==</em><b>true</b>, then
 *                <em>P<sub>j</sub></em>(<em>z</em>)<em> = pp</em>(<em>1,j</em>)<em>*S<sub>0</sub></em>(<em>z</em>)<em> + pp</em>(<em>2,j</em>)<em>*S<sub>1</sub></em>(<em>z</em>)<em> + ... + pp</em>(<em>n,j</em>)<em>*S<sub>n-1</sub></em>(<em>z</em>), where \n
 *                <em>S<sub>i</sub></em>(<em>z</em>) is the `translated' Chebyshev polynomial <em>S<sub>i</sub></em>((<em>z-c</em>)/<em>h</em>)<em> = T<sub>i</sub></em>(<em>z</em>), with
 *                <em>c = </em>(<em>intv</em>(<em>j</em>)<em>+intv</em>(<em>j+1</em>))<em> / 2</em> and
 *                <em>h = </em>(<em>intv</em>(<em>j+1</em>)<em>-intv</em>(<em>j</em>))<em> / 2</em>.
 *          </ul>
 *  \param  z0        is a number at which the piecewise polynomial <em>P</em>(<em>z</em>) is evaluated.
 *
 *  \return  The evaluated value of <em>P</em>(<em>z</em>) at <em>z=z0</em>.
 *
 *  \remark  If <em>z0</em> falls below the first interval, then the polynomial in the first interval will be used for evaluting <em>P</em>(<em>z0</em>).
 *  \remark  If <em>z0</em> flies over  the last  interval, then the polynomial in the last  interval will be used for evaluting <em>P</em>(<em>z0</em>).
 */
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
// note that if z0 falls below the first interval, then the polynomial in the first interval will be used for evaluting P(z0)
//           if z0 flies over  the last  interval, then the polynomial in the last  interval will be used for evaluting P(z0)
Real PiecewisePolynomialEvaluationInChebyshevBasis(const Matrix &pp, const Vector &intv, Real z0, bool basisTranslated=true);

//! Compute the weighted inner product of two piecewise polynomials expanded in a basis of `translated' (i.e. scale-and-shift) Chebyshev polynomials for each interval.
/*! \param  intervalWeights  defines the interval weights; <em>intervalWeights</em>(<em>j</em>) is the weight of the <em>j</em>th interval.
 *  \param  pp        is a matrix of Chebyshev coefficients which defines the piecewise polynomial <em>P</em>(<em>z</em>). \n
 *          For <em>z</em> in the <em>j</em>th interval, <em>P</em>(<em>z</em>) equals
 *          <em>P<sub>j</sub></em>(<em>z</em>) = <em>pp</em>(<em>1,j</em>)*<em>S<sub>0</sub></em>(<em>z</em>) + <em>pp</em>(<em>2,j</em>)*<em>S<sub>1</sub></em>(<em>z</em>)<em> + ... + pp</em>(<em>n,j</em>)<em>*S<sub>n-1</sub></em>(<em>z</em>).
 *  \param  qq        is a matrix of Chebyshev coefficients which defines the piecewise polynomial Q(z). \n
 *          For <em>z</em> in the <em>j</em>th interval, <em>Q</em>(<em>z</em>) equals
 *          <em>Q<sub>j</sub></em>(<em>z</em>) = <em>qq</em>(<em>1,j</em>)*<em>S<sub>0</sub></em>(<em>z</em>) + <em>qq</em>(<em>2,j</em>)*<em>S<sub>1</sub></em>(<em>z</em>)<em> + ... + qq</em>(<em>n,j</em>)<em>*S<sub>n-1</sub></em>(<em>z</em>).
 *
 *  \remark
 *  Here <em>S<sub>i</sub></em>(<em>z</em>) is the `translated' Chebyshev polynomial in that interval <em>[aa,bb]</em>,
 *     <em>S<sub>i</sub></em>((<em>z-c</em>)/<em>h</em>) = <em>T<sub>i</sub></em>(<em>z</em>),  <em>c</em> = (<em>aa+bb</em>))/<em>2</em>,  <em>h</em> = (<em>bb-aa</em>)/<em>2</em>,\n
 *  where <em>T<sub>i</sub></em>(<em>z</em>) is the Chebyshev polynomial of the first kind,
 *     <em>T<sub>0</sub></em>(<em>z</em>) = <em>1</em>, <em>T<sub>1</sub></em>(<em>z</em>) = <em>z</em>, and <em>T<sub>i</sub></em>(<em>z</em>) = <em>2*z*T<sub>i-1</sub></em>(<em>z</em>) - <em>T<sub>i-2</sub></em>(<em>Z</em>) for <em>i>=2</em>.
 *
 *  The (scaled) <em>j</em>th interval inner product is defined by
 *     <<em>P<sub>j</sub>,Q<sub>j</sub></em>> = (\f$\pi\f$/2) * [ <em>pp</em>(<em>1,j</em>)<em>*qq</em>(<em>1,j</em>) + sum<em><sub>k</sub> pp</em>(<em>k,j</em>)*<em>qq</em>(<em>k,j</em>) ],\n
 *  which comes from the property
 *     <<em>T<sub>0</sub>,T<sub>0</sub></em>>=\f$\pi\f$, <<em>T<sub>i</sub>,T<sub>i</sub></em>>=\f$\pi\f$/2 for <em>i>=1</em>, and <<em>T<sub>i</sub>,T<sub>j</sub></em>><em>=0</em> for <em>i\f$\neq\f$j</em>.
 *
 *  \return  The weighted inner product is <<em>P,Q</em>> = sum<em><sub>j</sub></em> [ <em>intervalWeights</em>(<em>j</em>) * <<em>P<sub>j</sub>,Q<sub>j</sub></em>> ].
 *
 *  \remark  For unit weights, pass an empty vector of <em>intervalWeights</em> (i.e. of length <em>0</em>).
 */
// compute the weighted inner product of two piecewise polynomials expanded
// in a basis of `translated' (i.e. scale-and-shift) Chebyshev polynomials for each interval
//
// pp and qq are two matrices of Chebyshev coefficients which define the piecewise polynomials P(z) and Q(z), respectively
// For z in the j-th interval, P(z) equals
//     P_j(z) = pp(1,j)*S_0(z) + pp(2,j)*S_1(z) + ... + pp(n,j)*S_{n-1}(z),
// and Q(z) equals
//     Q_j(z) = qq(1,j)*S_0(z) + qq(2,j)*S_1(z) + ... + qq(n,j)*S_{n-1}(z),
// where S_i(z) is the `translated' Chebyshev polynomial in that interval,
//     S_i((z-c)/h) = T_i(z),  c = (aa+bb)) / 2,  h = (bb-aa) / 2,
// with T_i(z) the Chebyshev polynomial of the first kind,
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
Real PiecewisePolynomialInnerProductInChebyshevBasis(const Matrix &pp, const Matrix &qq, const Vector &intervalWeights);

//! compute <em>Q</em>(<em>z</em>) = z*P(z), where <em>P</em>(<em>z</em>) and <em>Q</em>(<em>z</em>) are piecewise polynomials expanded in a basis of `translated' (i.e. scale-and-shift) Chebyshev polynomials for each interval.
/*! <em>P</em>(<em>z</em>) and <em>Q</em>(<em>z</em>) are stored as matrices of Chebyshev coefficients <em>pp</em> and <em>qq</em>, respectively.
 *  \param  intv      is a vector which defines the intervals. The <em>j</em>th interval is [<em>intv</em>(<em>j</em>), <em>intv</em>(<em>j+1</em>)).
 *  \param  pp        is a matrix of Chebyshev coefficients which defines the piecewise polynomial <em>P</em>(<em>z</em>). \n
 *          For <em>z</em> in the <em>j</em>th interval [<em>intv</em>(<em>j</em>), <em>intv</em>(<em>j+1</em>)), <em>P</em>(<em>z</em>) equals
 *          <em>P<sub>j</sub></em>(<em>z</em>) = <em>pp</em>(<em>1,j</em>)*<em>S<sub>0</sub></em>(<em>z</em>) + <em>pp</em>(<em>2,j</em>)*<em>S<sub>1</sub></em>(<em>z</em>)<em> + ... + pp</em>(<em>n,j</em>)<em>*S<sub>n-1</sub></em>(<em>z</em>).
 *  \return  A matrix <em>qq</em> of Chebyshev coefficients which defines the piecewise polynomial Q(z). \n
 *           For <em>z</em> in the <em>j</em>th interval [<em>intv</em>(<em>j</em>), <em>intv</em>(<em>j+1</em>)), <em>Q</em>(<em>z</em>) equals
 *           <em>Q<sub>j</sub></em>(<em>z</em>) = <em>qq</em>(<em>1,j</em>)*<em>S<sub>0</sub></em>(<em>z</em>) + <em>qq</em>(<em>2,j</em>)*<em>S<sub>1</sub></em>(<em>z</em>)<em> + ... + qq</em>(<em>n,j</em>)<em>*S<sub>n-1</sub></em>(<em>z</em>).
 *
 *  \remark
 *  Here <em>S<sub>i</sub></em>(<em>z</em>) is the `translated' Chebyshev polynomial in the <em>j</em>th interval, with
 *     <em>S<sub>i</sub></em>((<em>z-c</em>)/<em>h</em>) = <em>T<sub>i</sub></em>(<em>z</em>),  <em>c</em> = (<em>intv</em>(<em>j</em>) + <em>intv</em>(<em>j+1</em>))/<em>2</em>,  <em>h</em> = (<em><em>intv</em>(<em>j+1</em>) - <em>intv</em>(<em>j</em>)</em>)/<em>2</em>,\n
 *  where <em>T<sub>i</sub></em>(<em>z</em>) is the Chebyshev polynomial of the first kind,
 *     <em>T<sub>0</sub></em>(<em>z</em>) = <em>1</em>, <em>T<sub>1</sub></em>(<em>z</em>) = <em>z</em>, and <em>T<sub>i</sub></em>(<em>z</em>) = <em>2*z*T<sub>i-1</sub></em>(<em>z</em>) - <em>T<sub>i-2</sub></em>(<em>Z</em>) for <em>i>=2</em>.
 *
 *  \return  The matrix of coefficients <em>qq</em> which represents <em>Q</em>(<em>z</em>)<em> = z*P</em>(<em>z</em>).
 */
// compute Q(z) = z*P(z), where P(z) and Q(z) are piecewise polynomials expanded
// in a basis of `translated' (i.e. scale-and-shift) Chebyshev polynomials for each interval
//
// P(z) and Q(z) are stored as matrices of Chebyshev coefficients pp and qq, respectively
//
// for z in the j-th interval, P(z) equals
//     P_j(z) = pp(1,j)*S_0(z) + pp(2,j)*S_1(z) + ... + pp(n,j)*S_{n-1}(z),
// and Q(z) equals
//     Q_j(z) = qq(1,j)*S_0(z) + qq(2,j)*S_1(z) + ... + qq(n,j)*S_{n-1}(z),
// where S_i(z) is the `translated' Chebyshev polynomial in that interval,
//     S_i((z-c)/h) = T_i(z),  c = (intv(j)+intv(j+1))) / 2,  h = (intv(j+1)-intv(j)) / 2,
// with T_i(z) the Chebyshev polynomial of the first kind,
//     T_0(z) = 1, T_1(z) = z, and T_i(z) = 2*z*T_{i-1}(z) - T_{i-2}(z) for i>=2
//
// the returned matrix is qq which represents Q(z) = z*P(z)
Matrix PiecewisePolynomialInChebyshevBasisMultiplyX(const Matrix &pp, const Vector &intv);

/** @} */  // end of group_chebyshev



////////////////////////////////////////////////////////////////////////////
//    Conjugate Residual Method in the Polynomial Space
////////////////////////////////////////////////////////////////////////////

/** @defgroup group_cr_poly Conjugate Residual Method in Polynomial Space
 *
 *  @{
 */
//! This routine employs a conjugate-residual-type algorithm in polynomial space to minimize ||<em>P</em>(<em>z</em>)-<em>Q</em>(<em>z</em>)||<em><sub>w</sub></em>, where <em>P</em>(<em>z</em>), the base filter, is the input piecewise polynomial, and <em>Q</em>(<em>z</em>) is the output polynomial satisfying <em>Q</em>(<em>0</em>)<em>==1</em>, i.e. the constant term of <em>Q</em>(<em>z</em>) is <em>1</em>.
/*! Both <em>P</em>(<em>z</em>) and <em>Q</em>(<em>z</em>) are expanded in the `translated' (scale-and-shift) Chebyshev basis for each interval and
 *  presented as matrices of Chebyshev coefficients, denoted by <em>pp</em> and <em>qq</em>, respectively.
 *  \param  niter     is the number of conjugate-residual iterations. Therefore, the degree of <em>Q</em>(<em>z</em>) is up to <em>niter+1</em>.
 *  \param  intv      is a vector which defines the intervals.
 *          The <em>j</em>th interval is [<em>intv</em>(<em>j</em>),<em>intv</em>(<em>j+1</em>)).
 *  \param  w         is a vector of Chebyshev weights.
 *          The weight of <em>j</em>th interval is <em>w</em>(<em>j</em>). \n
 *          The interval weights define the inner product of two continuous functions and then
 *          the derived <em>w</em>-norm ||<em>P</em>(<em>z</em>)-<em>Q</em>(<em>z</em>)||<em><sub>w</sub></em>.
 *  \param  pp        is a matrix of Chebyshev coefficients which defines the piecewise polynomial <em>P</em>(<em>z</em>). \n
 *          For <em>z</em> in the <em>j</em>th interval [<em>intv</em>(<em>j</em>), <em>intv</em>(<em>j+1</em>)), <em>P</em>(<em>z</em>) equals
 *          <em>P<sub>j</sub></em>(<em>z</em>) = <em>pp</em>(<em>1,j</em>)*<em>S<sub>0</sub></em>(<em>z</em>) + <em>pp</em>(<em>2,j</em>)*<em>S<sub>1</sub></em>(<em>z</em>)<em> + ... + pp</em>(<em>niter+2,j</em>)<em>*S<sub>niter+1</sub></em>(<em>z</em>),
 *
 *  \remark  Here <em>S<sub>i</sub></em>(<em>z</em>) is the `translated' Chebyshev polynomial in the <em>j</em>th interval, with
 *           <em>S<sub>i</sub></em>((<em>z-c</em>)/<em>h</em>) = <em>T<sub>i</sub></em>(<em>z</em>),  <em>c</em> = (<em>intv</em>(<em>j</em>) + <em>intv</em>(<em>j+1</em>))/<em>2</em>,  <em>h</em> = (<em><em>intv</em>(<em>j+1</em>) - <em>intv</em>(<em>j</em>)</em>)/<em>2</em>,\n
 *           where <em>T<sub>i</sub></em>(<em>z</em>) is the Chebyshev polynomial of the first kind,
 *           <em>T<sub>0</sub></em>(<em>z</em>) = <em>1</em>, <em>T<sub>1</sub></em>(<em>z</em>) = <em>z</em>, and <em>T<sub>i</sub></em>(<em>z</em>) = <em>2*z*T<sub>i-1</sub></em>(<em>z</em>) - <em>T<sub>i-2</sub></em>(<em>Z</em>) for <em>i>=2</em>.
 *
 *  \return  A matrix, denoted by <em>qq</em>, represents a polynomial <em>Q</em>(<em>z</em>) with degree up to <em>1+niter</em>
 *           and satisfying <em>Q</em>(<em>0</em>)<em>==1</em>, such that ||<em>P</em>(<em>z</em>))-<em>Q</em>(<em>z</em>)||<em><sub>w</sub></em> is minimized. \n
 *           This polynomial <em>Q</em>(<em>z</em>) is expanded in the `translated' Chebyshev basis for each interval. \n
 *           To be precise, considering <em>z</em> in [<em>intv</em>(<em>j</em>)<em>,intv</em>(<em>j+1</em>)), <em>Q</em>(<em>z</em>) equals
 *           <em>Q<sub>j</sub></em>(<em>z</em>)<em> = qq</em>(<em>1,j</em>)<em>*S<sub>0</sub></em>(<em>z</em>)<em> + qq</em>(<em>2,j</em>)<em>*S<sub>1</sub></em>(<em>z</em>)<em> + ... + qq</em>(<em>niter+2,j</em>)<em>*S<sub>niter+1</sub></em>(<em>z</em>).
 *
 *  \remark  Since <em>Q</em>(<em>0</em>)<em>==1</em>, <em>P</em>(<em>0</em>)<em>==1</em> is expected.
 *           If <em>P</em>(<em>0</em>)\f$\neq\f$<em>1</em>, one can translate <em>P</em>(<em>z</em>).
 *           For example, if <em>P</em>(<em>0</em>)<em>==0</em>, one can use <em>1-P</em>(<em>z</em>) as input instead of <em>P</em>(<em>z</em>).
 *
 *  \remark  Typically, the base filter, defined by <em>pp</em> and <em>intv</em>, is from Hermite interpolation
 *           in intervals [<em>intv</em>(<em>j</em>)<em>,intv</em>(<em>j+1</em>)) for <em>j=1,...,nintv</em>,
 *           with <em>nintv</em> the number of intervals.
 */
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
Matrix FilteredConjugateResidualPolynomial(const Matrix &pp, const Vector &intv,
                                           const Vector &w, mkIndex niter);

//! This routine employs a conjugate-residual-type algorithm in polynomial space to compute <em>x=x0+s</em>(<em>A</em>)<em>*r0</em> with <em>r0=b-A*x0</em>, such that ||<em>1-z*s</em>(<em>z</em>)<em>-P</em>(<em>z</em>)||<em><sub>w</sub></em> is minimized.
/*! <ul>
 *  <li>  <em>P</em>(<em>z</em>) is a given piecewise polynomial, called the base filter. \n
 *        <em>P</em>(<em>z</em>) is expanded in the `translated' (scale-and-shift) Chebyshev basis for each interval,
 *        and presented as a matrix of Chebyshev coefficients <em>pp</em>.
 *  <li>  <em>s</em>(<em>z</em>) is a polynomial of degree up to niter, the number of conjugate-residual iterations.
 *  </ul>
 *  \param  A         is a sparse matrix.
 *  \param  x0,b      are vectors.
 *  \param  niter     is the number of conjugate-residual iterations. Therefore, the degree of <em>Q</em>(<em>z</em>) is up to <em>niter+1</em>.
 *  \param  intv      is a vector which defines the intervals.
 *          The <em>j</em>th interval is [<em>intv</em>(<em>j</em>),<em>intv</em>(<em>j+1</em>)).
 *  \param  w         is a vector of Chebyshev weights.
 *          The weight of <em>j</em>th interval is <em>w</em>(<em>j</em>). \n
 *          The interval weights define the inner product of two continuous functions and then
 *          the derived <em>w</em>-norm ||<em>P</em>(<em>z</em>)-<em>Q</em>(<em>z</em>)||<em><sub>w</sub></em>.
 *  \param  pp        is a matrix of Chebyshev coefficients which defines the piecewise polynomial <em>P</em>(<em>z</em>). \n
 *          For <em>z</em> in the <em>j</em>th interval [<em>intv</em>(<em>j</em>), <em>intv</em>(<em>j+1</em>)), <em>P</em>(<em>z</em>) equals
 *          <em>P<sub>j</sub></em>(<em>z</em>) = <em>pp</em>(<em>1,j</em>)*<em>S<sub>0</sub></em>(<em>z</em>) + <em>pp</em>(<em>2,j</em>)*<em>S<sub>1</sub></em>(<em>z</em>)<em> + ... + pp</em>(<em>niter+2,j</em>)<em>*S<sub>niter+1</sub></em>(<em>z</em>),
 *  \param  tol       is the tolerance; if the residual polynomial in <em>z</em>-norm is dropped by a factor lower
 *          than <em>tol</em>, then stop the conjugate-residual iteration.
 *
 *  \remark Here <em>S<sub>i</sub></em>(<em>z</em>) is the `translated' Chebyshev polynomial in the <em>j</em>th interval, with
 *         <em>S<sub>i</sub></em>((<em>z-c</em>)/<em>h</em>) = <em>T<sub>i</sub></em>(<em>z</em>),  <em>c</em> = (<em>intv</em>(<em>j</em>) + <em>intv</em>(<em>j+1</em>))/<em>2</em>,  <em>h</em> = (<em><em>intv</em>(<em>j+1</em>) - <em>intv</em>(<em>j</em>)</em>)/<em>2</em>,\n
 *         where <em>T<sub>i</sub></em>(<em>z</em>) is the Chebyshev polynomial of the first kind,
 *         <em>T<sub>0</sub></em>(<em>z</em>) = <em>1</em>, <em>T<sub>1</sub></em>(<em>z</em>) = <em>z</em>, and <em>T<sub>i</sub></em>(<em>z</em>) = <em>2*z*T<sub>i-1</sub></em>(<em>z</em>) - <em>T<sub>i-2</sub></em>(<em>Z</em>) for <em>i>=2</em>.
 *
 *  \return A vector <em>x=x0+s</em>(<em>A</em>)<em>*r0</em> with <em>r0=b-A*x0</em>,
 *         such that ||(<em>1-P</em>(<em>z</em>))<em>-z*s</em>(<em>z</em>)||<em><sub>w</sub></em> is minimized,
 *         subject to that <em>s</em>(<em>z</em>) is a polynomial of degree up to <em>niter</em>,
 *         where <em>P</em>(<em>z</em>) is the base filter.
 *         In short, <em>z*s</em>(<em>z</em>) approximates <em>1-P</em>(<em>z</em>).
 *
 *  \remark Since <em>z*s</em>(<em>z</em>) approximates <em>1-P</em>(<em>z</em>), <em>P</em>(<em>0</em>)<em>==1</em> is expected.
 *         If <em>P</em>(<em>0</em>)\f$\neq\f$<em>1</em>, one can translate <em>P</em>(<em>z</em>).
 *         For example, if <em>P</em>(<em>0</em>)<em>==0</em>, one can use <em>1-P</em>(<em>z</em>) as input instead of <em>P</em>(<em>z</em>).
 *
 *  \remark Typically, the base filter, defined by <em>pp</em> and <em>intv</em>, is from Hermite interpolation
 *         in intervals [<em>intv</em>(<em>j</em>)<em>,intv</em>(<em>j+1</em>)) for <em>j=1,...,nintv</em>,
 *         with <em>nintv</em> the number of intervals.
 *
 *  \remark A common application is to compute <em>R</em>(<em>A</em>)<em>*b</em>, where <em>R</em>(<em>z</em>) approximates <em>1-P</em>(<em>z</em>).
 *         In this case, one can set <em>x0 = 0</em> and then the return vector is <em>x = s</em>(<em>A</em>)<em>*b</em>, where
 *         <em>z*s</em>(<em>z</em>) approximates <em>1-P</em>(<em>z</em>); therefore, <em>A*x</em> is the wanted <em>R</em>(<em>A</em>)<em>*b</em>.
 *
 *  \see PolynomialFilterInterface::filteredSparseMatrixPolynomialVectorProduct().
 */
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
// see also PolynomialFilterInterface::filteredSparseMatrixPolynomialVectorProduct()
Vector FilteredConjugateResidualMatrixPolynomialVectorProduct(const SparseMatrix &A, const Vector &x0, const Vector &b,
                                           const Matrix &pp, const Vector &intv, const Vector &w,
                                           mkIndex niter, Real tol=0.0);

//! An interface to define <em>R</em>(<em>S</em>), i.e. set <em>R</em> and <em>S</em>, and compute <em>R</em>(<em>S</em>)<em>*v</em>, where <em>R</em> is the residual polynomial, <em>S</em> is a sparse matrix, and <em>v</em> is a vector.
/*! <ul>
 *  <li>  <em>S</em> is a sparse symmetric matrix.
 *  <li>  <em>R</em>(<em>z</em>) is a polynomial which approximates a piecewise polynomial <em>1-P</em>(<em>z</em>), with <em>P</em>(<em>z</em>) the base filter.
 *  <li>  <em>v</em> is the input vector.
 *  </ul>
 *
 *  In short, use setFilter() to define <em>R</em>(<em>S</em>) once, and
 *  use filteredSparseMatrixPolynomialVectorProduct() to compute <em>R</em>(<em>S</em>)<em>*v</em> for every <em>v</em>.
 */
// an interface to define R(S), i.e. set R and S, and compute R(S)*v, where
// S is a sparse symmetric matrix,
// R(z) is a polynomial which approximates a piecewise polynomial 1-P(z), with P(z) the base filter, and
// v is the input vector
//
// in short, use setFilter() to define R(S) once, and
// use filteredSparseMatrixPolynomialVectorProduct() to compute R(S)*v for every v
class PolynomialFilterInterface {
private:
    // works for SparseMatrix only
    //! A pointer to the input sparse symmetric matrix via setFilter().
    // a pointer to the input sparse symmetric matrix via setFilter()
    static const SparseMatrix *Sptr;
    //! The sparse symmetric matrix S in computing <em>R</em>(<em>S</em>)<em>*v</em>; <em>S</em> is translated from the input sparse symmetric matrix <em>S0</em> in setFilter().
    // the sparse symmetric matrix S in computing R(S)*v; S is translated from the input sparse symmetric matrix S0 in setFilter()
    static SparseMatrix S;

    // parameters for polynomial filtering
    //! The <em>j</em>th interval is [<em>intervals</em>(<em>j</em>),<em>intervals</em>(<em>j+1</em>)).
    /*! The intervals are decided by GetIntervals() invoked by setFilter() which takes parameters <em>frame</em>, <em>baseDeg</em>, and <em>polyDeg</em>:
     *  <ul>
     *  <li>  <em>frame</em> is a vector consisting of 4 ordered elements:
     *        [<em>frame</em>(<em>1</em>)<em>,frame</em>(<em>4</em>)] is the interval which (tightly) contains all eigenvalues of <em>S</em>, and
     *        [<em>frame</em>(<em>2</em>)<em>,frame</em>(<em>3</em>)] is the interval in which the eigenvalues are sought.
     *  <li>  <em>baseDeg</em> is the left-and-right degree of the base filter for each interval.
     *  <li>  <em>polyDeg</em> is the (maximum possible) degree of <em>s</em>(<em>z</em>), with <em>z*s</em>(<em>z</em>) the polynomial filter.
     *  </ul>
     */
    // the j-th interval is [intervals(j),intervals(j+1))
    // the intervals are decided by GetIntervals() invoked by setFilter() which takes parameters frame, baseDeg, and polyDeg:
    // 1. frame is a vector consisting of 4 ordered elements:
    //     [frame(1),frame(4)] is the interval which (tightly) contains all eigenvalues of S, and
    //     [frame(2),frame(3)] is the interval in which the eigenvalues are sought
    // 2. baseDeg is the left-and-right degree of the base filter for each interval
    // 3. polyDeg is the (maximum possible) degree of z*s(z), with s(z) the polynomial filter
    static Vector intervals;
    //! <em>intervalWeights</em>(<em>j</em>) is the weight of <em>j</em>th interval [<em>intervals</em>(<em>j</em>),<em>intervals</em>(<em>j+1</em>)).
    // intervalWeights(j) is the weight of j-th interval [intervals(j), intervals(j+1))
    static Vector intervalWeights;
    //! Maximum possible degree of <em>s</em>(<em>z</em>), with <em>z*s</em>(<em>z</em>) the polynomial filter.
    // maximum possible degree of s(z), with z*s(z) the polynomial filter
    static mkIndex polyDegree;
    //! Left-and-right degree of the base filter in each interval.
    // left-and-right degree of the base filter in each interval
    static mkIndex baseDegree;
    //! A base filter, which is a piecewise polynomial typically from Hermite interpolation.
    /*! For <em>j</em>th interval, the polynomial is expanded in the `translated' (scale-and-shift) Chebyshev basis,
     *  with the Chebyshev coefficients stored in <em>baseFilter</em>(:,<em>j</em>).
     */
    // a base filter, which is a piecewise polynomial typically from Hermite interpolation
    // for j-th interval, the polynomial is expanded in the `translated' (scale-and-shift) Chebyshev basis,
    // with the Chebyshev coefficients stored in baseFilter(:,j)
    static Matrix baseFilter;
    //! A zero vector of length the number of rows/columns of <em>S</em>.
    // a zero vector of length the number of rows/columns of S
    static Vector zn;

public:
    //! The information of computing the polynomial filter. See class PolynomialFilterInfo in polyfilt.h for more information.
    // the information of computing the polynomial filter; see class PolynomialFilterInfo in polyfilt.h for more information
    static PolynomialFilterInfo filterInfo;

    // routines for the polynomial filtered Lanczos eigensolver
    //! This member function sets the sparse matrix <em>S</em> translated from <em>S0</em>, and also the base filter <em>P</em>(<em>z</em>).
    /*! \param  S0 is the sparse symmetric matrix of concern.
     *  \param  frame is a vector of <em>4</em> ordered elements. \n
     *          [<em>frame</em>(<em>1</em>),<em>frame</em>(<em>4</em>)] is the interval which (tightly) contains all eigenvalues of <em>S0</em>, and
     *          [<em>frame</em>(<em>2</em>),<em>frame</em>(<em>3</em>)] is the interval in which the eigenvalues of <em>S0</em> are sought.
     *  \param  polyDeg is the (maximum possible) degree of s(z), with z*s(z) the polynomial filter.
     *  \param  baseDeg is the left-and-right degree of the base filter for each interval.
     *  \param  opts is a collection of options to determine the intervals.
     *
     *  \see GetIntervals().
     */
    // this member function sets the sparse matrix S translated from S0, and also the base filter P(z)
    // S0 is the sparse symmetric matrix of concern
    // frame is a vector of 4 ordered elements
    //     [frame(1),frame(4)] is the interval which (tightly) contains all eigenvalues of S0, and
    //     [frame(2),frame(3)] is the interval in which the eigenvalues of S0 are sought
    // polyDeg is the (maximum possible) degree of s(z), with z*s(z) the polynomial filter
    // baseDeg is the left-and-right degree of the base filter for each interval
    // opts is a collection of options to determine the intervals
    static Vector setFilter(const SparseMatrix &S0, const Vector &frame, mkIndex polyDeg, mkIndex baseDeg, IntervalOptions &opts);
    //! This routine computes <em>R</em>(<em>S</em>)<em>*v</em>, where it is assumed that <em>S</em> and <em>R</em>(<em>S</em>) have been defined by setFilter().
    /*! \remark  <em>R</em>(<em>z</em>) is in the form <em>z*s</em>(<em>z</em>) which minimizes ||(<em>1-P</em>(<em>z</em>))-<em>z*s</em>(<em>z</em>)||<em><sub>w</sub></em> by a conjugate-residual-type algorithm,
     *           where <em>P</em>(<em>z</em>) is a piecewise polynomial as the base filter, and <em>s</em>(<em>z</em>) is a polynomial of degree up to PolynomialFilterInterface::polyDegree,
     *           and <em>w</em>-norm is defined by data members PolynomialFilterInterface::intervals and PolynomialFilterInterface::intervalWeights.
     *  \remark  <em>P</em>(<em>z</em>) is defined by the data members PolynomialFilterInterface::baseFilter and PolynomialFilterInterface::intervals.
     *           It is expanded in the `translated' (scale-and-shift) Chebyshev basis for each interval, and presented as a matrix of Chebyshev coefficients PolynomialFilterInterface::baseFilter.
     *
     *  \see FilteredConjugateResidualMatrixPolynomialVectorProduct().
     */
    // this routine computes R(S)*v, where it is assumed that S and R(z) have been defined by setFilter()
    // R(z) is in the form z*s(z) which minimizes ||(1-P(z))-z*s(z)||_w by a conjugate-residual-type algorithm,
    // where P(z) is a piecewise polynomial as the base filter, and s(z) is a polynomial of degree up to polyDegree,
    // and w-norm is defined by data members intervals and intervalWeights
    //
    // note that P(z) is defined by the data members baseFilter and intervals
    // it is expanded in the `translated' (scale-and-shift) Chebyshev basis for each interval,
    // and presented as a matrix of Chebyshev coefficients baseFilter
    //
    // see also FilteredConjugateResidualMatrixPolynomialVectorProduct()
    static Vector filteredSparseMatrixPolynomialVectorProduct(const Vector &v) {
        return (*Sptr)*FilteredConjugateResidualMatrixPolynomialVectorProduct(*Sptr, zn, v, baseFilter, intervals, intervalWeights, polyDegree);
    }
};

/** @} */  // end of group_cr_poly


#endif
