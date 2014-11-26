#ifndef BASICS_H
#define BASICS_H
/*! \file basics.h
 *  The following constants and static template functions are provided:
 *  <ul>
 *  <li> The machine epsilon \f$\epsilon_M\f$ and the machine infinity \f$\infty\f$.
 *  <li> The mathematical constants \f$\pi\f$ and the Euler number \f$e\f$.
 *  <li> Square \f$x^2\f$, absolute value \f$|x|\f$, minimum \f$\min(x,y)\f$ and maximum \f$\max(x,y)\f$.
 *  </ul>
 */

#include <stdlib.h>  // for exit, atof, atoi, atol, srand, rand, etc.
#include <iostream>  // for ostream, etc.


#ifdef USE_NAMESPACE
namespace MATKIT {
#endif


/*! \brief Class Basics defines the machine constants, and mathematical constants in numerical form,
 *  and also providing a few basic static template functions.
 */
// class Basics defines the machine constants, and mathematical constants in numerical form, as well as providing a few basic static template functions
class Basics {
public:
    // some basics constants
    //! The difference between 1 and the smallest exactly representable number greater than 1 in floating-point arithmetic.
    /*! In IEEE 754 floating-point arithmetic,
     *  the machine epsilon is <em>2<sup>-52</sup></em> in double precision format,
     *  or <em>2<sup>-23</sup></em> in single precision format.
     */
    // the difference between 1 and the smallest exactly representable number greater than 1 in floating-point arithmetic
    // in IEEE 754 floating-point arithmetic, the machine epsilon is 2^{-52} in double precision format, or 2^{-23} in single precision format
    static const Real machineEpsilon;
    //! The largest possible value in the floating-point representation.
    /*! In IEEE 754 floating-point arithmetic,
     *  the largest possible value is (<em>2-2<sup>-52</sup></em>)*<em>2<sup>1023</sup></em> in double precision format,
     *  or (<em>2-2<sup>-23</sup></em>)*<em>2<sup>127</sup></em> in single precision format.
     */
    // the largest possible value in the floating-point representation
    // in IEEE 754 floating-point arithmetic, the largest possible value is (2-2^{-52})*2^{1023} in double precision format, or (2-2^{-23})*2^{127} in single precision format
    static const Real infinity;
    //! The mathematical constant 3.14159265..., the ratio of any circle's circumference to its diameter.
    // the mathematical constant 3.14159265..., the ratio of any circle's circumference to its diameter
    static const Real pi;
    //! The mathematical constant 2.718281828..., often referred to as the Euler number.
    // the mathematical constant 2.718281828..., often referred to as the Euler number
    static const Real naturalLog;

    // output ostream for error messages
    //! A pointer to an ostream for the output error messages. By default, it points to <b>std::cerr</b>.
    // a pointer to an ostream for the output error messages; by default, it points to std::cerr
    static std::ostream *err;

    // square and abs
    //! Return <em>x<sup>2</sup></em>, the square of <em>x</em>.
    // return x^2, the square of x
    template<class Element>
    static Element square(Element x) { return x*x; }
    //! Return |<em>x</em>|, the absolute value of <em>x</em>.
    // return |x|, the absolute value of x
    template<class Element>
    static Element abs(Element x) {
        return (x>=0.0 ? x: -x);
    }

    // max and min
    //! Return the larger value of <em>x</em> and <em>y</em>.
    // return the larger value of x and y
    template<class OrderedElement>
    static OrderedElement max(OrderedElement x, OrderedElement y) {
        return (x>=y ? x: y);
    }
    //! Return the smaller value of <em>x</em> and <em>y</em>.
    // return the smaller value of x and y
    template<class OrderedElement>
    static OrderedElement min(OrderedElement x, OrderedElement y) {
        return (x<=y ? x: y);
    }

    //! Quit the program.
    // quit the program
    static void quit(int status) {
        #ifdef USE_MEX
        mexErrMsgTxt("Quit the mex-function due to an error!");
        #else
        exit(status);
        #endif
    }
};


#ifdef USE_NAMESPACE
}  // end of namespace MATKIT
#endif

#endif  // end of #ifndef BASICS_H
