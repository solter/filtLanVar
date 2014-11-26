#ifndef MATKITDEF_H
#define MATKITDEF_H
/*! \file matkitdef.h
 *  Three fundamental parameters are decided here: <b>Real</b>, <tt>USE_BLAS</tt>/<tt>USE_MWBLAS</tt>, and <tt>USE_MEX</tt>.
 *  <ul>
 *  <li>  <b>Real</b> is normally either <b>double</b> (by <tt>typedef double Real;</tt>) or
 *        <b>single</b> (by <tt>typedef single Real;</tt>).
 *
 *  <li>  All <b>MATKIT</b> routines are implemented, and the <b>BLAS</b> routines can be <em>optionally</em> invoked
 *        by setting USE_BLAS to improve the computational efficiency.
 *        For example, if <b>Real</b> is <b>double</b> and <tt>USE_BLAS</tt> is set
 *        by either `<tt>\#define USE_BLAS</tt>' or by the compiler flag `<tt>-DUSE_BLAS</tt>', then
 *        the <b>BLAS</b> routines <b>dscale</b>_() <b>dcopy</b>_(), <b>daxpy</b>_(), <b>dnrm2</b>_(), <b>ddot</b>_(), <b>dger</b>_(),
 *        <b>dspr</b>_(), <b>dspr2</b>_(), <b>dgemv</b>_(), <b>dspmv</b>_(), <b>dgemm</b>_()
 *        are declared and invoked in various <b>MATKIT</b> routines.
 *
 *  <li>  The <b>MATKIT</b> routines can be interfaced to <b>OCTAVE</b>/<b>MATLAB</b>, by setting <tt>USE_MEX</tt>.
 *  </ul>
 *  See source code for details.
 */

/* The following parameters are now handled by the compiler parameters "-DUSE_NAMESPACE", "-DUSE_SINGLE", "-DUSE_BLAS", "-DUSE_SPBLAS", and "-DUSE_MEX", "-DUSE_MWBLAS"
   #define USE_NAMESPACE
   #define USE_SINGLE
   #define USE_BLAS
   #define USE_SPBLAS
   #define USE_MEX
   #define USE_MWBLAS
*/

#ifdef USE_MEX
    // for mex routines
    #include "mex.h"

    // for defining an ostream "mexostream" for printing messages to OCTAVE/MATLAB console
    #include <stdio.h>    // for EOF
    #include <iostream>   // for ostream, streambuf
#endif


#ifdef USE_NAMESPACE
namespace MATKIT {
#endif


//! Define <b>Real</b> numbers, normally either <b>double</b> or <b>single</b>.
/*! One may set <tt>USE_SINGLE</tt> by either `<tt>\#define USE_SINGLE</tt>' in matkitdef.h or the compiler flag `<tt>-DUSE_SINGLE</tt>'
 *  to define <b>Real</b> as <b>single</b>.
 */
#ifdef USE_SINGLE
    typedef float Real;
#else
    typedef double Real;
#endif

//! Define <b>mkIndex</b> for index values.
/*! <b>mkIndex</b> is typically <b>int</b>, <b>unsigned</b>, <b>long</b>, or <b>unsigned long</b>.
 *  However, to invoke <b>Intel MKL SPBLAS</b>, <b>mkIndex</b> must be <b>int</b> or <b>unsigned</b>.
 */
// mkIndex is typically int, unsigned, long, or unsigned long
// to invoke Intel MKL SPBLAS, mkIndex must be int or unsigned
typedef unsigned mkIndex;

//! Define <b>mkSignedIndex</b> for signed index values.
/*! <b>mkSignedIndex</b> must be signed, typically <b>int</b> or <b>long</b>.
 */
// mkSignedIndex must be signed, typically int or long
typedef int mkSignedIndex;


#if defined(USE_MWBLAS)
    #include "blas.h"
#elif defined(USE_BLAS)
    extern "C" {
        // BLAS routines
        #ifdef USE_SINGLE
        // in-place scale (or a matrix) x by alp, i.e. x = alp*x
        void sscal_(const int *len, const float *alp, const float *x, const int *incx);
        // copy a vector (or a matrix) x to a vector (or a matrix) y
        void scopy_(const int *len, const float *x, const int *incx, float *y, const int *incy);
        // compute y = alp*x + y, where alp is a constant, and x, y are vectors (or matrices)
        void saxpy_(const int *len, const float *alp, const float *x, const int *incx, const float *y, const int *incy);
        // 2-norm of a vector, or Frobenius norm of a matrix
        float snrm2_(const int *len, const float *x, const int *incx);
        // inner (dot) product x'*y
        float sdot_(const int *len, const float *x, const int *incx, const float *y, const int *incy);
        // outer product A = alp*x*y'
        void sger_(const int *m, const int *n, const Real *alp, const float *x, const int *incx, const float *y, const int *incy, float *aa, const int *lda);
        // outer product A = alp*x*x', where A is stored as a symmetric matrix (in packed form)
        void sspr_(const char *uplo, const int *n, const float *alp, const float *x, const int *incx, float *ap);
        // symmetric outer product A = alp*(x*y'+y*x'), where A is stored as a symmetric matrix (in packed form)
        void sspr2_(const char *uplo, const int *n, const float *alp, const float *x, const int *incx, const float *y, const int *incy, float *ap);
        // matrix - vector product
        void sgemv_(const char *trans, const int *m, const int *n, const float  *alp, const float  *aa, const int *lda, const float  *x, const int *incx, const float  *beta, float  *y, const int *incy);
        // symmetric matrix - vector product
        void sspmv_(const char *uplo,                const int *n, const float  *alp, const float  *ap,                 const float  *x, const int *incx, const float  *beta, float  *y, const int *incy);
        // matrix - matrix multiplication
        void sgemm_(const char *transA, const char *transB, const int *m, const int *n, const int *k, const float  *alp, const float  *aa, const int *lda, const float  *bb, const int *ldb, const float  *beta, float  *cc, const int *ldc);
        #else
        // in-place scale (or a matrix) x by alp, i.e. x = alp*x
        void dscal_(const int *len, const double *alp, const double *x, const int *incx);
        // copy a vector (or a matrix) x to a vector (or a matrix) y
        void dcopy_(const int *len, const double *x, const int *incx, double *y, const int *incy);
        // compute y = alp*x + y, where alp is a constant, and x, y are vectors (or matrices)
        void daxpy_(const int *len, const double *alp, const double *x, const int *incx, const double *y, const int *incy);
        // 2-norm of a vector, or Frobenius norm of a matrix
        double dnrm2_(const int *len, const double *x, const int *incx);
        // inner (dot) product
        double ddot_(const int *len, const double *x, const int *incx, const double *y, const int *incy);
        // outer product A = alp*x*y'
        void dger_(const int *m, const int *n, const Real *alp, const double *x, const int *incx, const double *y, const int *incy, double *aa, const int *lda);
        // symmetric outer product A = alp*(x*y'+y*x'), where A is stored as a symmetric matrix (in packed form)
        void dspr2_(const char *uplo, const int *n, const double *alp, const double *x, const int *incx, const double *y, const int *incy, double *ap);
        // outer product A = alp*x*x', where A is stored as a symmetric matrix (in packed form)
        void dspr_(const char *uplo, const int *n, const double *alp, const double *x, const int *incx, double *ap);
        // matrix - vector product
        void dgemv_(const char *trans, const int *m, const int *n, const double *alp, const double *aa, const int *lda, const double *x, const int *incx, const double *beta, double *y, const int *incy);
        // symmetric matrix - vector product
        void dspmv_(const char *uplo,                const int *n, const double *alp, const double *ap,                 const double *x, const int *incx, const double *beta, double *y, const int *incy);
        // matrix - matrix multiplication
        void dgemm_(const char *transA, const char *transB, const int *m, const int *n, const int *k, const double *alp, const double *aa, const int *lda, const double *bb, const int *ldb, const double *beta, double *cc, const int *ldc);
        #endif
    }
#endif

#ifdef USE_SPBLAS
    extern "C" {
        // Intel MKL SPBLAS routines
        // compute y = alpha*A*x + beta*y
        #ifdef USE_SINGLE
        void mkl_scscmv_(const char *transa, int *nr, int *nc, float *alpha,
                         const char *matdescra, const float *sa, const int *rowIdx, const int *colPtr, const int *colPtr1,
                         const float *v, float *beta, float *w);
        void mkl_scsrmv_(const char *transa, int *nr, int *nc, float *alpha,
                         const char *matdescra, const float *sa, const int *rowIdx, const int *colPtr, const int *colPtr1,
                         const float *v, float *beta, float *w);
        #else
        void mkl_dcscmv_(const char *transa, int *nr, int *nc, double *alpha,
                         const char *matdescra, const double *sa, const int *rowIdx, const int *colPtr, const int *colPtr1,
                         const double *v, double *beta, double *w);
        void mkl_dcsrmv_(const char *transa, int *nr, int *nc, double *alpha,
                         const char *matdescra, const double *sa, const int *rowIdx, const int *colPtr, const int *colPtr1,
                         const double *v, double *beta, double *w);
        #endif
    }
#endif


#ifdef USE_MEX
// define an ostream "mexostream" for printing messages to MATLAB console
// References:
// [1] Prototypes of std::streambuf and std::ostream:
// http://www.cplusplus.com/reference/iostream/streambuf/
// http://www.cplusplus.com/reference/iostream/ostream/
// [2] A working code (subject to minor changes):
// http://stackoverflow.com/questions/243696/correctly-over-loading-a-stringbuf-to-replace-cout-in-a-matlab-mex-file
class mexbuf: public std::streambuf {
public:
    mexbuf() {}

protected:
    virtual int_type overflow(int_type ch = EOF) {
        // in a test, this function is called whenever "<< std::endl" and the input ch is '\n'
        if (ch != EOF) {
            mexPrintf("%.1s", &ch);
        }
        return 1;
    }

    virtual std::streamsize xsputn(const char *str, std::streamsize num) {
        mexPrintf("%.*s", num, str);
        return num;
    }
};

// class mexostream
class mexostream: public std::ostream {
protected:
    mexbuf buf;

public:
    mexostream(): std::ostream(&buf) {}
};

// declare mexcout as an instance of mexostream
extern mexostream mexcout;
// mexcout is defined in matkitfunc.cpp
#endif


#ifdef USE_NAMESPACE
}  // end of namespace MATKIT
#endif

#endif  // end of #ifndef MATKITDEF_H
