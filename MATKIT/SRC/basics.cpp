//============================================================================
// the MATKIT class Basics
// last update August 2011
//============================================================================

#include <string.h>  // for strlen, memcpy, etc.
#include <math.h>    // for HUGE_VAL

#include "matkitdef.h"
#include "basics.h"

#ifdef USE_NAMESPACE
namespace MATKIT {
#endif


// use IEEE 754 floating-point arithmetic
#ifdef USE_SINGLE
    const Real Basics::machineEpsilon = 1.19209289e-7;           // 2^(-23)
    const Real Basics::infinity =       3.402824e+38;            // (2-2^(-23))*2^127
#else
    const Real Basics::machineEpsilon = 2.220446049250313e-16;   // 2^(-52)
    // const Real Basics::infinity =    1.797693134862316e+308;  // (2-2^(-52))*2^1023
    const Real Basics::infinity =       HUGE_VAL;
#endif

const Real Basics::pi = 3.1415926535897931;
const Real Basics::naturalLog = 2.7182818284590455;


#ifdef USE_MEX
std::ostream *Basics::err = &mexcout;
#else
std::ostream *Basics::err = &std::cerr;
#endif


#ifdef USE_NAMESPACE
}  // end of namespace MATKIT
#endif
