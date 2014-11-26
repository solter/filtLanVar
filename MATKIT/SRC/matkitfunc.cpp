//============================================================================
// the MATKIT routines w/o classes
// coded by H.-r. Fang
// last update April, 2012
//============================================================================

#include <stdlib.h>  // for exit, atof, atoi, atol, srand, rand, etc.
#include <string.h>  // for memcpy, strcmp, strncmp, etc.
#include <stddef.h>  // for NULL pointer, pointer subtraction, etc.
#include <math.h>    // for sqrt, pow, etc.
#include <iostream>  // for cout, cerr, endl, etc. (under namespace std)

#include <stdio.h>   // for fopen, fscanf, fclose, etc.
#include <ctype.h>   // for tolower, toupper, etc.

#include "matkitdef.h"
#include "matkitfunc.h"

using std::endl;

#ifdef USE_MEX
// for overriding operators "new" and "delete"
#include <new>
#endif


#ifdef USE_NAMESPACE
namespace MATKIT {
#endif

#ifdef USE_MEX
// for printing messages to MATLAB console
mexostream mexcout;
#endif


////////////////////////////////////////////////////////////////////////////
//    memory routines
////////////////////////////////////////////////////////////////////////////

#ifdef USE_MEX_MALLOC
// the overridden new and delete operators invokes mxMalloc() and mxFree()
// however, in tests under OCTAVE ver 3.2.4 and 3.4.2, and MATLAB ver 2010b and 2011a,
// the resulting code with the overridden new and delete operators is unstable (possibly due to a double free)
// at this point it is suggested not to define "USE_MEX_MALLOC"
// the drawback without overriding the new and delete operators is that,
// when mexErrMsgTxt() is called (e.g. in Basics::quit()), memory leaks are expected
// http://www.devx.com/tips/Tip/5608

// override operators "new"
void *operator new(size_t sz) throw (std::bad_alloc) {
    void *ptr = mxMalloc((mwSize)sz);
    if (ptr == NULL)
        throw std::bad_alloc();
    return ptr;
}
void *operator new(size_t sz, const std::nothrow_t &nothrow_constant) throw() {
    return mxMalloc((mwSize)sz);
}
void *operator new[](size_t sz) throw (std::bad_alloc) {
    void *ptr = mxMalloc((mwSize)sz);
    if (ptr == NULL)
        throw std::bad_alloc();
    return ptr;
}
void *operator new[](size_t sz, const std::nothrow_t &nothrow_constant) throw() {
    return mxMalloc((mwSize)sz);
}

// override operators "delete"
void operator delete(void *p) throw() {
    if (p != NULL)
        mxFree(p);
    // found no documentation about whether it is safe to mxFree a null pointer
    // note however in a test under MATLAB ver 2010b (code compiled by mex),
    //    deleting a null pointer is handled without getting into this routine;
    //    hence the line "if (p != NULL)" above should be unnecessary
    // no clue whether it is compiler dependent
    // so "if (p != NULL)" would be better added
    // the comment applies to the following 3 routines
}
void operator delete(void *p, const std::nothrow_t &nothrow_constant) throw() {
    if (p != NULL)
        mxFree(p);
}
void operator delete[](void *p) throw() {
    if (p != NULL)
        mxFree(p);
}
void operator delete[](void *p, const std::nothrow_t &nothrow_constant) throw() {
    if (p != NULL)
        mxFree(p);
}

// References:
// [1] C++ The Core Language (chapter 7) by Gregory Satir and Doug Brown, published by O'Reilly Media, 1995.
// [2] Prototypes of C++ new/delete operators:
// http://www.cplusplus.com/reference/std/new/operator%20new/
// http://www.cplusplus.com/reference/std/new/operator%20new%5B%5D/
// http://www.cplusplus.com/reference/std/new/operator%20delete/
// http://www.cplusplus.com/reference/std/new/operator%20delete%5B%5D/
// [3] A working code by Petter Strandmark 2009:
// http://www.mathworks.com/matlabcentral/fileexchange/24349-use-c-new-and-delete-in-matlab-mex-files&watching=24349
#endif



////////////////////////////////////////////////////////////////////////////
//    miscellaneous routines
////////////////////////////////////////////////////////////////////////////

// compute des0[j*inc2] = src0[j*inc1] for j=0,...,len-1, where
//    src0=src if inc1>=0, or src0=src-(len-1)*inc1 if inc1<0, and
//    des0=des if inc2>=0, or des0=des-(len-1)*inc2 if inc2<0
// conceptually, this routine computes w=v, where
//    v is the input  vector stored in src0[0,inc1,...,(len-1)*inc1], and
//    w is the output vector stored in des0[0,inc2,...,(len-1)*inc2]
// it can be used to compute A=B, where A and B are both general or both symmetric matrices
// it can be used to copy nonzero elements of a sparse matrix
// it can also be used to compute a submatrix or a subvector (i.e. copy elements)
// this routine invokes the level 1 BLAS routine xCOPY if USE_BLAS is defined
void memory_xcopy(mkIndex len, const Real *src, Real *des, mkSignedIndex inc1, mkSignedIndex inc2) {
// MWBLAS routines scopy() and dcopy() require non-const *src, so it cannot be used in memory_xcopy(),
// which asks const *src
#ifdef USE_BLAS
    int len2 = (int)len;
    int incx = (int)inc1;
    int incy = (int)inc2;
    #ifdef USE_SINGLE
        scopy_(&len2, src, &incx, des, &incy);
    #else
        dcopy_(&len2, src, &incx, des, &incy);
    #endif
#else
    if ((inc1==1 && inc2==1) || (inc1==-1 && inc2==-1)) {
        memcpy(des, src, len*sizeof(Real));
    }
    else {
        if (inc1 < 0) {
            src -= ((mkSignedIndex)len*inc1-inc1);
        }
        if (inc2 < 0) {
            des -= ((mkSignedIndex)len*inc2-inc2);
        }
        while (len--) {
            *des = *src;
            src += inc1;
            des += inc2;
        }
    }
#endif
}

// compute sv[i*inc] *= alp for i=0,...,len-1
// conceptually, this routine computes v *= alp, where v is stored in sv[0,inc,...,(len-1)*inc]
// it can also be used to compute A *= alp, where A is a general, symmetric, or sparse matrix
// this routine invokes the level 1 BLAS routine xSCAL if USE_BLAS is defined
void in_place_scale_elements(mkIndex len, Real alp, Real *sv, mkIndex inc) {
#if defined(USE_MWBLAS)
    mwSignedIndex len2 = (mwSignedIndex)len;
    mwSignedIndex inc2 = (mwSignedIndex)inc;
    #ifdef USE_SINGLE
        sscal(&len2, &alp, sv, &inc2);
    #else
        dscal(&len2, &alp, sv, &inc2);
    #endif
#elif defined(USE_BLAS)
    // use to scale the elements
    int len2 = (int)len;
    int inc2 = (int)inc;
    #ifdef USE_SINGLE
        sscal_(&len2, &alp, sv, &inc2);
    #else
        dscal_(&len2, &alp, sv, &inc2);
    #endif
#else
    // scale the elements
    if (alp != 1.0) {
        while (len--) {
            (*sv) *= alp;
            sv += inc;
        }
    }
#endif
}

// compute su[i] = sv[i] + bet*sw[i] (if su!=sv) or sw[i] += bet*sv[i] (if su==sv) for i=0,...,len-1
// conceptually, this routine computes u = v+bet*w (if su!=sv) or v += bet*w (if su==sv), where
//    v, w are the input vectors of length len stored in sv[], wv[], respectively, and
//    u is the output vector of length len stored in su[]
// if su==NULL is given, then the required memory of size len will be allocated for the output vector u;
//    otherwise it is assumed that sufficient memory has been allocated, with the address pointed to by su
// this routine can be used for computing C = A+bet*B or A += bet*B, where
//    A and B are both general or both symmetric matrices
// this routine invokes the level 1 BLAS routine xAXPY if USE_BLAS is defined
void elementwise_addition(mkIndex len, const Real *sv, Real bet, const Real *sw, Real *&su) {
    if (len == 0) {
        // an exception, empty input vectors, do nothing
        return;
    }

    if (su == NULL) {
        // allocate memory
        su = new Real[len];
    }

    if (su != sv) {
        // compute u = v
        memory_xcopy(len, sv, su);
    }

    if (bet != 0.0) {
        // compute u = u + bet*w
        #ifdef USE_BLAS
            int len2 = (int)len;
            int inc = 1;
            #ifdef USE_SINGLE
                saxpy_(&len2, &bet, sw, &inc, su, &inc);
            #else
                daxpy_(&len2, &bet, sw, &inc, su, &inc);
            #endif
        #else
            mkIndex z = len;
            Real *u0 = su;
            const Real *w0 = sw;
            // compute su[i] += alp*sw[i] for i=0,...,len-1
            while (z--)
                *(u0++) += bet*(*w0++);
        #endif
    }
}


#ifdef USE_MEX
// convert an input argument mxArray of OCTAVE/MATLAB to a double precision number
// always output "double" even if "Real" is "float", to reduce loss of significant digits due to data-type conversion
double mxArray_to_double(const mxArray *ma, mkIndex offset) {
    if (!mxIsNumeric(ma))
        mexErrMsgTxt("mxArray2SingleNumber(const mxArray *): the input is not numeric");
    void *val = mxGetData(ma);

    switch (mxGetClassID(ma)) {
        case mxDOUBLE_CLASS:
            return (Real)(*((double *)(val)+offset));
        case mxSINGLE_CLASS:
            return (Real)(*((float *)(val)+offset));
        case mxINT8_CLASS:
            return (Real)(*((int8_t *)(val)+offset));
        case mxUINT8_CLASS:
            // return (Real)(*((uint8_t *)(val)+offset));  // "uint8_t" not recognized by OCTAVE compiler mkoctfile ver 3.2.4
            return (Real)(*((unsigned char *)(val)+offset));
        case mxINT16_CLASS:
            return (Real)(*((int16_t *)(val)+offset));
        case mxUINT16_CLASS:
            // return (Real)(*((uint16_t *)(val)+offset));  // "uint16_t" not recognized by OCTAVE compiler mkoctfile ver 3.2.4
            return (Real)(*((unsigned short *)(val)+offset));
        case mxINT32_CLASS:
            return (Real)(*((int32_t *)(val)+offset));
        case mxUINT32_CLASS:
            // return (Real)(*((uint32_t *)(val)+offset));  // "uint32_t" not recognized by OCTAVE compiler mkoctfile ver 3.2.4
            return (Real)(*((unsigned *)(val)+offset));
        default:
            // the remains are not numeric (cell, string, or structure) and has been handled in "if (!mxIsNumeric(ma)" above
            // should not end in here unless new types of mxArray are introduced
            return 0.0;
    }
}
#endif  // of #ifdef USE_MEX



////////////////////////////////////////////////////////////////////////////
//    vector routines
////////////////////////////////////////////////////////////////////////////

// print out vector v to ostream s, where v is stored in sv[] of length len
void print_vector(std::ostream &s, mkIndex len, const Real *sv) {
    std::ostream *sptr = &s;
    #ifdef USE_MEX
    if (&s == &std::cout || &s == &std::cerr) {
        // stdout and stderr are not supported by the mex interface
        // use mexostream (which inherits ostream) "mexcout" for printing messages to OCTAVE/MATLAB console
        // see matkitdef.h for more information
        sptr = &mexcout;
    }
    #endif

    *sptr << "    ";  // the leading spaces
    if (len == 0) {
        // a special case, empty vector
        *sptr << "  []";
    }
    for (mkIndex i=0; i<len; i++) {
        // print v(i+1)
        *sptr << "  " << sv[i];
    }
    *sptr << endl;
}

// write a vector v to a plain text file fileName[],
//    where v is of length len stored in sv[]
// the format of the file:
//    each line contains a number
//    the first line has the length of the vector as an integer
//    after that, the (i+1)st has has the i-th element of the matrix as a Real number, for i=1,...,len
// the return value:
//    0: a successful write
//    1: file open error
//    2: file close error
//    a negative integer: a number returned by fprintf() which signals an error
int vector_write(mkIndex len, const Real *sv, const char fileName[]) {
    // open the file to write
    FILE *mmFile = fopen(fileName, "w");
    if (mmFile == NULL) {
        // file open error!
        return 1;
    }

    // print the length
    int info;
    info = fprintf(mmFile, "%u\n", (unsigned)len);
    if (info < 0)
        return info;

    // print the elements
    for (mkIndex i=0; i<len; i++) {
        #ifdef USE_SINGLE
        info = fprintf(mmFile, "%.7g\n",   sv[i]);
        #else
        info = fprintf(mmFile, "%.15lg\n", sv[i]);
        #endif
        if (info < 0) {
            // file write error
            fclose(mmFile);
            return info;
        }
    }

    // now close the file
    if (fclose(mmFile) == EOF) {
        // close file error
        return 2;
    }

    // a successful write, return 0
    return 0;
}

// read a vector from a plain text file fileName[] with the format:
//    each line contains a number
//    the first line has the length of the vector as an integer
//    after that, the (i+1)st has has the i-th element of the matrix as a Real number, for i=1,...,len
// the elements of the vector are stored in sv[] and the length is recorded as len
//    
// the return value:
//    0: a successful read
//   -1: file open error
//   -2: file is empty
//   -3: length of the vector error
//   -4: the vector data is incomplete or invalid
int vector_read(const char fileName[], mkIndex &len, Real *&sv) {
    // open the file to read
    FILE *mmFile = fopen(fileName, "r");
    if (mmFile == NULL) {
        // file open error!
        return -1;
    }

    char string[256];  // assume that there are at most 255 characters in each line of the header
    if (fgets(string, 256, mmFile) == NULL) {
        // file is empty!
        fclose(mmFile);
        return -2;
    }

    // get rid of the comment/empty lines
    while (strlen(string)==0 || *string=='\n' || *string=='%') {
        if (fgets(string, 256, mmFile) == NULL) {
            // no information about the length of the vector!?
            fclose(mmFile);
            return -3;
        }
    }

    // read the length of the vector
    int len0 = atoi(string);
    if (len0 < 0) {
        // length is negative!?
        fclose(mmFile);
        return -3;
    }
    len = (mkIndex)len0;

    // allocate memory if necessary
    if (len && sv == NULL) {
        sv = new Real[len];
    }

    // read the vector
    Real *v0 = sv;
    while (len0--) {
        #ifdef USE_SINGLE
        if (fscanf(mmFile, "%g\n",  v0++) != 1) {
        #else
        if (fscanf(mmFile, "%lg\n", v0++) != 1) {
        #endif
            // incomplete data!?
            fclose(mmFile);
            return -4;
        }
    }

    // a successful read, close the input file and return 0
    fclose(mmFile);
    return 0;
}

// compute w = v {op} alp (if sw!=sv), or v {op}= alp (if sw==sv), where
// {op} == +,-,*,/ if flag == 0,1,2,3, respectively, where
//    alp is a scalar
//    v is the input  vector of length len stored in sv[]
//    w is the output vector of length len stored in sw[]
// if sw==NULL is given, then the required memory of size len will be allocated for the output vector w;
//    otherwise, it is assumed that sufficient memory has been allocated, with the address pointed to by sw
// this routine invokes in_place_scale_elements() if sw == sv and flag == 2 (i.e. v *= alp)
// this routine can also be used to compute B = A{op}alp (with sv == sw) and A {op}= alp (with sv != sw), where
//    A and B are both general or both symmetric matrices
// this routine does not catch the division-by-zero exception (i.e. alp == 0.0 and flag == 3);
//    a division-by-zero will return an NaN (not-a-number) for each element
void vector_scalar_operation(mkIndex len, const Real *sv, Real alp, Real *&sw, unsigned flag) {
    if (len == 0) {
        // an exception, empty input, do nothing
        return;
    }

    if (sw == NULL) {
        // allocate memory
        sw = new Real[len];
    }

    Real *w0 = sw;
    if (sw == sv) {
        if (flag == 0) {
            // w[i] += alp for i=0,...,len-1
            if (alp != 0.0) {
                while (len--)
                    *w0++ += alp;
            }
        }
        else if (flag == 1) {
            // w[i] -= alp for i=0,...,len-1
            if (alp != 0.0) {
                while (len--)
                    *w0++ -= alp;
            }
        }
        else if (flag == 2) {
            // w[i] *= alp for i=0,...,len-1
            in_place_scale_elements(len, alp, sw);
        }
        else if (flag == 3) {
            // w[i] /= alp for i=0,...,len-1
            if (alp != 1.0) {
                while (len--)
                    *w0++ /= alp;
                // one could call "in_place_scale_elements(len, 1.0/alp, sw)" instead
            }
        }
    }
    else {
        if (flag == 0) {
            // w[i] = v[i] + alp for i=0,...,len-1
            if (alp != 0.0) {
                while (len--)
                    *(w0++) = *(sv++) + alp;
            }
            else
                memcpy(sw, sv, len*sizeof(Real));
        }
        else if (flag == 1) {
            // w[i] = v[i] - alp for i=0,...,len-1
            if (alp != 0.0) {
                while (len--)
                    *(w0++) = *(sv++) - alp;
            }
            else
                memcpy(sw, sv, len*sizeof(Real));
        }
        else if (flag == 2) {
            // w[i] = v[i] * alp for i=0,...,len-1
            if (alp != 1.0) {
                while (len--)
                    *(w0++) = *(sv++) * alp;
            }
            else
                memcpy(sw, sv, len*sizeof(Real));

        }
        else if (flag == 3) {
            // w[i] = v[i] * alp for i=0,...,len-1
            if (alp != 1.0) {
                while (len--)
                    *(w0++) = *(sv++) / alp;
            }
            else
                memcpy(sw, sv, len*sizeof(Real));
        }
    }
}

// form a (len+|k|)-by-(len+|k|) general matrix A with the k-th diagonal formed by vector v and otherwise zero
// v is the input vector of length len stored in sv[]
// A is the output (len+|k|)-by-(len+|k|) general matrix A stored in sa[], such that
//    if k >= 0, then A(j,j+k)=v(j) for j=1,...,len and otherwise zero, or
//    if k <  0, then A(j-k,j)=v(j) for j=1,...,len and otherwise zero
// if sa==NULL is given, then the required memory of size (len+|k|)*(len+|k|) will be allocated for the output matrix A;
//    otherwise it is assumed that sufficient memory has been allocated, with the address pointed to by sa
// this routine resembles A = diag(v,k) in OCTAVE/MATLAB
void vector_to_diag(mkIndex len, const Real *sv, Real *&sa, mkSignedIndex k) {
    // dimension of the output (square) matrix
    mkIndex nrc = (k>=0) ? (len+k) : (len-k);

    if (nrc == 0) {
        // an exception, empty output, do nothing
        return;
    }

    // memory size of the output matrix
    mkIndex sz = nrc*nrc;

    if (sa == NULL) {
        // allocate memory
        sa = new Real[sz];
    }

    // zero initialization
    Real *a0 = sa;
    while (sz--)
        *a0++ = 0.0;

    mkIndex shift = nrc + 1u;
    if (k >= 0) {
        // for a pointer to A(1,1+k)
        a0 = sa + k*nrc;
    }
    else {
        // for a pointer to A(1-k,1)
        a0 = sa - k;
    }

    // copy sv[0,...,len-1] to the k-th diagonal of A
    while (len--) {
        *a0 = *sv++;
        a0 += shift;
    }
}

// form a (len+|k|)-by-(len+|k|) sparse matrix A with the k-th diagonal formed by a vector v and otherwise zero
// v is the input vector of length len stored in sv[]
// A is the output (len+|k|)-by-(len+|k|) sparse matrix A stored in store[], rowIdx[], colPtr[], such that
//    if k >= 0, A(j,j+k)=v(j) for j=1,...,len and otherwise zero, or
//    if k <  0, A(j-k,j)=v(j) for j=1,...,len and otherwise zero
// if store==NULL (rowIdx==NULL, or colPtr==NULL) is given, then the required memory will be allocated;
//    otherwise it is assumed that sufficient memory has been allocated, with the address pointed to by store (rowIdx, or colPtr, respectively)
// this routine resembles A = spdiag(v,k) in OCTAVE, a sparse version of diag(v,k) in OCTAVE/MATLAB
void vector_to_spdiag(mkIndex len, const Real *sv, Real *&sa, mkIndex *&rowIdx, mkIndex *&colPtr, mkSignedIndex k) {
    // dimension of the output (square) matrix
    mkIndex nrc = (k>=0) ? (len+k) : (len-k);

    // find the number of nonzero elements
    mkIndex nnz = 0;
    for (mkIndex i=0; i<len; i++) {
        if (sv[i] != 0.0)
            nnz ++;
    }

    // allocate memory
    if (sa == NULL && nnz)
        sa = new Real[nnz];
    if (rowIdx == NULL && nnz)
        rowIdx = new mkIndex[nnz];
    if (colPtr == NULL)
        colPtr = new mkIndex[nrc+1u];

    // set the diagonal matrix
    mkIndex z = 0;
    if (k >= 0) {
        // the diagonal is in the upper triangular part
        // set A(:,1:k) to be all zero
        for (mkIndex i=0; i<(mkIndex)k; i++)
             colPtr[i] = 0;

        // set A(:,k+1:nrc) be with A(j,j+k) = sv[j-1] for j=1,...,len and otherwise zero, where nrc = len+k
        for (mkIndex i=0; i<len; i++) {
            colPtr[i+k] = z;
            if (sv[i] != 0.0) {
                rowIdx[z] = i;
                sa[z++] = sv[i];
            }
        }
        // now z == nnz
        colPtr[nrc] = nnz;
    }
    else {  // k < 0
        // the diagonal is in the strictly lower triangular part
        // set A(:,1:len) be with A(j-k,j) = sv[j-1] for j=1,...,len and otherwise zero
        for (mkIndex i=0; i<len; i++) {
            colPtr[i] = z;
            if (sv[i] != 0.0) {
                rowIdx[z] = i-k;
                sa[z++] = sv[i];
            }
        }
        // now z == nnz

        // set A(:,len+1:nrc) to be all zero, where nrc = len-k
        for (mkIndex i=len; i<=nrc; i++)
            colPtr[i] = nnz;
    }
}


#ifdef USE_MEX
// convert *mxArray to a vector
void mxArray_to_vector(const mxArray *mv, mkIndex &len, Real *&sv) {
    if (!mxIsNumeric(mv))
        mexErrMsgTxt("mxArray_to_vector(const mxArray*, mkIndex&, Real*&): the input mxArray* contains non-numeric data!");
    if (mxIsSparse(mv))
        mexErrMsgTxt("mxArray_to_vector(const mxArray*, mkIndex&, Real*&): the input mxArray* contains sparse data!");

    // OCTAVE/MATLAB: [nr, nc] = size(v);
    mkIndex nr = mxGetM(mv);
    mkIndex nc = mxGetN(mv);

    // exception handling
    if (nr > 1 && nc > 1)
        mexErrMsgTxt("mxArray_to_vector(const mxArray*, mkIndex&, Real*&): the input *mxArray has multiple rows and multiple columns, not a vector!");

    if (nr == 0 || nc == 0) {
        // in case of no element, length of the vector should be zero
        len = 0;
        return;
    }

    // OCTAVE/MATLAB: len = length(v);
    len = (nr>=nc) ? nr : nc;

    if (sv == NULL) {
        // allocate memory
        sv = new Real[len];
    }

    // copy the data
    if (mxIsDouble(mv)) {
        double *vv = mxGetPr(mv);
        #ifdef USE_SINGLE
            mexWarnMsgTxt("mxArray_to_vector(const mxArray*, mkIndex&, Real*&): the input mxArray* contains double precision floating-point data, but the output will be in single precision!");
            mkIndex z = len;
            Real *v0 = sv;
            while (z--)
                *v0++ = (Real)(*vv++);
        #else
            memory_xcopy(len, vv, sv);
        #endif
    }
    else if (mxIsSingle(mv)) {
        float *vv = (float *)mxGetData(mv);
        #ifndef USE_SINGLE
            mexWarnMsgTxt("mxArray_to_vector(const mxArray*, mkIndex&, Real*&): the input mxArray* contains single precision floating-point data, but the output will be in double precision!");
            mkIndex z = len;
            Real *v0 = sv;
            while (z--)
                *v0++ = (Real)(*vv++);
        #else
            memory_xcopy(len, vv, sv);
        #endif
    }
    else {
        mexErrMsgTxt("mxArray_to_vector(const mxArray*, mkIndex&, Real*&): the input mxArray* contains integer-type data; cannot convert it!");
    }
    if (mxIsComplex(mv))
        mexWarnMsgTxt("mxArray_to_vector(const mxArray*, mkIndex&, Real*&): the input mxArray* is a complex vector; the imaginary part is ignored!");
}

// convert a vector (of length len) to *mxArray (of size len-by-1)
mxArray *vector_to_mxArray(mkIndex len, const Real *sv) {
    // create the mxArray
    mxArray *mv = mxCreateDoubleMatrix(len, 1, mxREAL);

    // copy the data
    double *store = mxGetPr(mv);
    #ifdef USE_SINGLE
        const Real *v0 = sv;
        while (len--)
            *store++ = (double)(*v0++);
    #else
        memory_xcopy(len, sv, store);
    #endif

    return mv;
}
#endif  // of #ifdef USE_MEX



////////////////////////////////////////////////////////////////////////////
//    matrix routines
////////////////////////////////////////////////////////////////////////////

//! Print out matrix <em>A</em> to ostream <em>s</em>, where <em>A</em> is of size <em>nr</em>-by-<em>nc</em> stored in <em>sa</em>[].
// print out matrix A to ostream s, where A is of size nr-by-nc stored in sa[]
void print_matrix(std::ostream &s, mkIndex nr, mkIndex nc, const Real *sa) {
    std::ostream *sptr = &s;
    #ifdef USE_MEX
    if (&s == &std::cout || &s == &std::cerr) {
        // stdout and stderr are not supported by the mex interface
        // use mexostream (which inherits ostream) "mexcout" for printing messages to OCTAVE/MATLAB console
        // see matkitdef.h for more information
        sptr = &mexcout;
    }
    #endif

    if (nr == 0 && nc == 0) {
        // a special case, empty matrix
        *sptr << "      []" << endl;
    }
    for (mkIndex i=0; i<nr; i++) {
        // print A(i+1,:)
        *sptr << "    ";
        for (mkIndex j=0; j<nc; j++) {
             // print A(i+1,j+1)
             *sptr << "  " << sa[nr*j+i];
        }
        *sptr << endl;
    }
}

// convert a general matrix A to a sparse matrix B, such that mathematically A == B, where
//    A is the nr-by-nc input general matrix stored in sa[], and
//    B is the nr-by-nc output sparse matrix in CSC format stored in sb[], rowIdx[], colPtr[]
// if patternOnly==true (default is false), then only the sparsity pattern is considered, and sb[] will not be used
// if sb==NULL (rowIdx==NULL, or colPtr==NULL) is given, then the required memory will be allocated;
//    otherwise it is assumed that sufficient memory has been allocated, with the address pointed to by sb (rowIdx, or colPtr, respectively)
// this routine resembles B = sparse(A) in OCTAVE/MATLAB
void full_to_sparse(mkIndex nr, mkIndex nc, const Real *sa, Real *&sb, mkIndex *&rowIdx, mkIndex *&colPtr, bool patternOnly) {
    if (colPtr == NULL) {
        // allocate memory for column pointers
        colPtr = new mkIndex[nc+1u];
    }

    // for pointer to A(1,1)
    const Real *a0 = sa;

    if ((sb == NULL && !patternOnly) || rowIdx == NULL) {
        // find the number of nonzero elements
        mkIndex nnz = 0;
        mkIndex z = nr*nc;
        while (z--) {
            if (*a0++ != 0.0)
                nnz ++;
        }
        if (nnz) {
            // allocate memory if necessary
            if (sb == NULL && !patternOnly)
                sb = new Real[nnz];
            if (rowIdx == NULL)
                rowIdx = new mkIndex[nnz];
        }
        // reset the pointer to A(1,1)
        a0 = sa;
    }

    // now convert full A to sparse B
    Real *b0 = sb;
    mkIndex *ridx = rowIdx;
    mkIndex *cptr = colPtr;
    // set colPtr[0] = 0
    *cptr = 0;
    for (mkIndex j=1; j<=nc; j++) {
        // counter for nonzero elements in A(:,j)
        mkIndex m = 0;

        // OCTAVE/MATLAB: B(:,j) = A(:,j)
        for(mkIndex i=0; i<nr; i++) {
            if (*a0 != 0.0) {
                if (!patternOnly)
                    *b0++ = *a0;
                *ridx++ = i;
                m ++;
            }
            a0 ++;
        }
        // now m is the number of nonzero elements in A(:,j)

        // set colPtr[j] = colPtr[j-1] + m
        *(cptr+1) = *cptr + m;

        // for pointer to colPtr[j]
        cptr ++;
    }
}

// compute a vector v formed by the k-th diagonal of a general matrix A
// to be precise, this routine computes v(i)=A(i,i+k) if k>=0, or v(i)=A(i-k,i) if k<0, where
//    A is the nr-by-nc input general matrix stored in sa[]
//    v is the output vector v stored in sv[] containing the k-th diagonal elements (in order) of A
// the return value is the length of v
// if sv==NULL is given, then the required memory will be allocated for the output vector v;
//    otherwise, it is assumed that sufficient memory has been allocated, with the address pointed to by sv
// this routine resembles v = diag(A,k) in OCTAVE/MATLAB
mkIndex matrix_diag(mkIndex nr, mkIndex nc, const Real *sa, Real *&sv, mkSignedIndex k) {
    // compute the length of output vector and
    // a pointer to the first element of A to copy
    mkIndex len;
    const Real *a0;
    if (k >= 0) {
        len = (nc>(mkIndex)k) ? (nc-(mkIndex)k) : 0;
        if (len > nr)
            len = nr;
        if (len == 0)
            return 0;
        // pointer to A(1,k)
        a0 = sa + k*nr;
    }
    else {  // k < 0
        len = (nr>(mkIndex)(-k)) ? (nr-(mkIndex)(-k)) : 0;
        if (len > nc)
            len = nc;
        if (len == 0)
            return 0;
        // pointer to A(1-k,1)
        a0 = sa - k;
    }

    if (sv == NULL) {
        // allocate memory
        sv = new Real[len];
    }

    // copy elements
    mkIndex shift = nr + 1u;
    for (mkIndex i=0; i<len; i++) {
        sv[i] = *a0;
        a0 += shift;
    }

    // return the length of output vector
    return len;
}

// transpose a general matrix A and store the result as a general matrix B (not in-place)
// to be precise, mathematically B(j,i)=A(i,j) for i=1,...,nr and j=1,...,nc, where
//    A is the nr-by-nc input matrix A stored in sa[]
//    B, the transpose of A, is the nc-by-nr output matrix B stored in sb[]
// if sb==NULL is given, then the required memory of size nr*nc will be allocated for the output general matrix B;
//    otherwise, it is assumed that sufficient memory has been allocated, with the address pointed to by sb
// this routine resembles B = A' or equivalently B = transpose(A) in OCTAVE/MATLAB
void matrix_transpose(mkIndex nr, mkIndex nc, const Real *sa, Real *&sb) {
    if (nr == 0 || nc == 0) {
        // a special case, empty input
        return;
    }

    if (sb == NULL) {
        // allocate memory
        sb = new Real[nr*nc];
    }

    if (nr > nc) {
        // addresses of A(1,1) and B(1,1)
        const Real *a1 = sa;
        Real *b0 = sb;

        // transpose the rows of A to be columns of B
        mkIndex ii = nr;
        while (ii--) {
            // compute B(:,i) = A(i,:), where i = nr-ii
            const Real *a0 = a1++;
            // a0 is the address of A(i,1), b0 is the address of B(1,i), and
            // a1 is the address of A(i+1,1) now
            mkIndex jj = nc;
            while (jj--) {
                *b0++ = *a0;
                a0 += nr;
            }
        }
    }
    else {
        // addresses of A(1,1) and B(1,1)
        const Real *a0 = sa;
        Real *b1 = sb;

        // transpose the columns of A to be rows of B
        mkIndex jj = nc;
        while (jj--) {
            // compute B(j,:) = A(:,j), where j = nc-jj
            Real *b0 = b1++;
            // a0 is the address of A(j,1), b0 is the address of B(1,j), and
            // b1 is the address of B(1,j+1) now
            mkIndex ii = nr;
            while (ii--) {
                *b0 = *a0++;
                b0 += nc;
            }
        }
    }
}

// compute a general matrix B as the submatrix formed by elements in rows i1,...,i2 and columns j1,...,j2 of a general matrix A
// 1. the input is an nr-by-nc general matrix A stored in sa[]
// 2. the output is B=A(i1:i2,j1:j2) if indexFrom1==true, or B=A(i1+1:i2+1,j1+1:j2+1) if indexFrom1==false,
//    where B is a general matrix stored in sb[]
// the row indices i1,i2 and column indices j1,j2 range over 1,...,nr and 1,...,nc (if indexFrom1==true), or
//    over 0,...,nr-1 and 0,...,nc-1 (if indexFrom1==false), respectively
// remarks:
// 1. if sb==NULL is given, then the required memory of size nr*nc will be allocated for the output matrix B;
//    otherwise it is assumed that sufficient memory has been allocated, with the address pointed to by sb
// 2. this routine allows i1>i2 and/or j1>j2
//    if i1>i2 and j2>j2, then the output is B=A(i1:-1:i2,j1:-1:j2) if indexFrom1==true, or B=A(i1+1:-1:i2+1,j1+1:-1:j2+1) if indexFrom1==false
void submatrix_of_matrix(mkIndex nrows, mkIndex ncols, const Real *sa,
                         mkIndex i1, mkIndex i2, mkIndex j1, mkIndex j2, Real *&sb, bool indexFrom1) {
    // find dimensions of the submatrix
    mkIndex nr, nc;
    mkSignedIndex stride;
    if (i1 <= i2) {
        nr = i2 - i1 + 1u;
        if (j1 <= j2) {
            nc = j2 - j1 + 1u;
            stride = nrows - nr;
        }
        else {  // j1 > j2
            nc = j1 - j2 + 1u;
            stride = -(mkSignedIndex)nrows - (mkSignedIndex)nr;
        }
    }
    else {  // i1 > i2
        nr = i1 - i2 + 1u;
        if (j1 <= j2) {
            nc = j2 - j1 + 1u;
            stride = nrows + nr;
        }
        else {  // j1 > j2
            nc = j1 - j2 + 1u;
            stride = -(mkSignedIndex)nrows + (mkSignedIndex)nr;
        }
    }

    if (sb == NULL) {
        // allocate memory
        sb = new Real[nr*nc];
    }

    Real *b0 = sb;
    const Real *a0 = sa;
    if (indexFrom1)
        a0 += ((i1-1) + (j1-1)*nrows);
    else
        a0 += (i1 + j1*nrows);
    // now b0 is the address of B(1,1), and a0 is the address of A(i1,j1)
    // we assume indexFrom1==true in the comments hereafter in this routine

    // copy elements
    mkIndex jj = nc;
    if (i1 <= i2) {
        while (jj--) {
            // compute B(:,j+1) = A(i1:i2,j1+s*j), where j = nc-jj-1, and
            // s=1 if j1<=j2, or s=-1 if j1>j2
            mkIndex ii = nr;
            while (ii--)
                *b0++ = *a0++;
            // for the address of A(i1,j1+s*(j+1))
            a0 += stride;
        }
    }
    else {  // i1 > i2
        while (jj--) {
            // compute B(:,j+1) = A(i1:-1:i2,j1+s*j), where j = nc-jj-1, and
            // s=1 if j1<=j2, or s=-1 if j1>j2
            mkIndex ii = nr;
            while (ii--)
                *b0++ = *a0--;
            // for the address of A(i1,j1+s*(j+1))
            a0 += stride;
        }
    }
}

// permute columns of a general matrix A (not in-place)
// the result is a general matrix B, such that
//    B(i,perm[j-1])   = A(i,j) for i=1,...,nr and j=1,...,nc (if indexFrom1==true), or
//    B(i,perm[j-1]+1) = A(i,j) for i=1,...,nr and j=1,...,nc (if indexFrom1==false), where
//    A and B are the nr-by-nc input and output matrices stored in sa[] and sb[], respectively
// the indices perm[0,...,nc-1] range over 1,...,nc (if indexFrom1==true), or over 0,...,nc-1 (if indexFrom1==false)
// if sb==NULL is given, then the required memory of size nr*nc will be allocated for the output matrix B;
//    otherwise it is assumed that sufficient memory has been allocated, with the address pointed to by sb
// this routine resembles B(:,perm) = A in OCTAVE/MATLAB
void permute_matrix_columns(mkIndex nr, mkIndex nc, const mkIndex *perm, const Real *sa, Real *&sb, bool indexFrom1) {
    if (nr == 0 || nc == 0) {
        // a special case, empty input, do nothing
        return;
    }

    if (sb == NULL) {
        // allocate memory
        sb = new Real[nr*nc];
    }

    Real *b0 = sb;
    if (indexFrom1) {
        // pre-shift the pointer
        b0 -= nr;
    }
    for (mkIndex j=0; j<nc; j++) {
        mkIndex j2 = perm[j];
        // compute B(:,j2+1) = A(:,j+1) if indexFrom1==false,
        //      or B(:,j2)   = A(:,j+1) if indexFrom1==true
        memory_xcopy(nr, sa, b0+j2*nr);
        sa += nr;
    }
}

// permute rows of a general matrix A (not in-place)
// the result is a general matrix B, such that
//    B(perm[i-1],j)   = A(i,j) for i=1,...,nr and j=1,...,nc (if indexFrom1==true), or
//    B(perm[i-1]+1,j) = A(i,j) for i=1,...,nr and j=1,...,nc (if indexFrom1==false), where
//    A and B are the nr-by-nc input and output matrices stored in sa[] and sb[], respectively
// the indices perm[0,...,nr-1] range over 1,...,nr (if indexFrom1==true), or over 0,...,nr-1 (if indexFrom1==false)
// if sb==NULL is given, then the required memory of size nr*nc will be allocated for the output matrix B;
//    otherwise it is assumed that sufficient memory has been allocated, with the address pointed to by sb
// this routine resembles B(perm,:) = A in OCTAVE/MATLAB
void permute_matrix_rows(mkIndex nr, mkIndex nc, const mkIndex *perm, const Real *sa, Real *&sb, bool indexFrom1) {
    if (nr == 0 || nc == 0) {
        // a special case, empty input, do nothing
        return;
    }

    if (sb == NULL) {
        // allocate memory
        sb = new Real[nr*nc];
    }

    Real *b0 = sb;
    if (indexFrom1) {
        // pre-shift the pointer
        b0 --;
    }
    for (mkIndex i=0; i<nr; i++) {
        mkIndex i2 = perm[i];
        // compute B(i2+1,:) = A(i+1,:) if indexFrom1==false,
        //      or B(i2  ,:) = A(i+1,:) if indexFrom1==true
        memory_xcopy(nc, sa, b0+i2, nr, nr);
        sa ++;
    }
}

#ifdef USE_MEX
// convert *mxArray to a matrix
void mxArray_to_matrix(const mxArray *mA, mkIndex &nr, mkIndex &nc, Real *&sa) {
    if (!mxIsNumeric(mA))
        mexErrMsgTxt("mxArray_to_matrix(const mxArray*, mkIndex&, mkIndex&, Real*&): the input mxArray* contains non-numeric data!");
    if (mxIsSparse(mA))
        mexErrMsgTxt("mxArray_to_matrix(const mxArray*, mkIndex&, mkIndex&, Real*&): the input mxArray* contains sparse data!");

    // OCTAVE/MATLAB: [nr, nc] = size(A);
    nr = mxGetM(mA);
    nc = mxGetN(mA);

    if (nr == 0 || nc == 0) {
        // a special case, empty case, do nothing
        return;
    }

    // copy the data
    mkIndex sz = nr*nc;
    if (sa == NULL)
        sa = new Real[sz];
    if (mxIsDouble(mA)) {
        double *vv = mxGetPr(mA);
        #ifdef USE_SINGLE
            mexWarnMsgTxt("mxArray_to_matrix(const mxArray*, mkIndex&, mkIndex&, Real*&): the input mxArray* contains double precision floating-point data, but the output will be in single precision!");
            mkIndex z = sz;
            Real *a0 = sa;
            while (z--)
                *a0++ = (Real)(*vv++);
        #else
            memory_xcopy(sz, vv, sa);
        #endif
    }
    else if (mxIsSingle(mA)) {
        float *vv = (float *)mxGetData(mA);
        #ifndef USE_SINGLE
            mexWarnMsgTxt("mxArray_to_matrix(const mxArray*, mkIndex&, mkIndex&, Real*&): the input mxArray* contains single precision floating-point data, but the output will be in double precision!");
            mkIndex z = sz;
            Real *a0 = sa;
            while (z--)
                *a0++ = (Real)(*vv++);
        #else
            memory_xcopy(sz, vv, sa);
        #endif
    }
    else {
        mexErrMsgTxt("mxArray_to_matrix(const mxArray*, mkIndex&, mkIndex&, Real*&): the input mxArray* contains integer-type data; cannot convert it!");
    }
    if (mxIsComplex(mA))
        mexWarnMsgTxt("mxArray_to_matrix(const mxArray*, mkIndex&, mkIndex&, Real*&): the input mxArray* is a complex matrix; the imaginary part is ignored!");
}

// convert a matrix to *mxArray
mxArray *matrix_to_mxArray(mkIndex nr, mkIndex nc, const Real *sa) {
    // create the mxArray
    mxArray *mA = mxCreateDoubleMatrix(nr, nc, mxREAL);

    // copy the data
    double *store = mxGetPr(mA);
    mkIndex sz = nr*nc;
    #ifdef USE_SINGLE
        while (sz--)
            *store++ = (double)(*sa++);
    #else
        memory_xcopy(sz, sa, store);
    #endif

    return mA;
}
#endif  // of #ifdef USE_MEX



////////////////////////////////////////////////////////////////////////////
//    symmetric matrix routines
////////////////////////////////////////////////////////////////////////////

// print out symmetric matrix A to ostream s, where A is of size nrc-by-nrc stored in packed form in sa[]
void print_symmatrix(std::ostream &s, mkIndex nrc, const Real *sa) {
    std::ostream *sptr = &s;
    #ifdef USE_MEX
    if (&s == &std::cout || &s == &std::cerr) {
        // stdout and stderr are not supported by the mex interface
        // use mexostream (which inherits ostream) "mexcout" for printing messages to OCTAVE/MATLAB console
        // see matkitdef.h for more information
        sptr = &mexcout;
    }
    #endif

    if (nrc == 0) {
        // a special case, empty matrix
        *sptr << "      []" << endl;
    }
    for (mkIndex i=0; i<nrc; i++) {
        s << "    ";
        for (mkIndex j=0; j<nrc; j++) {
            // print A(i+1,j+1)
            if (i >= j)
                s << "  " << sa[j*(2u*nrc-j+1u)/2u+(i-j)];
            else
                s << "  " << sa[i*(2u*nrc-i+1u)/2u+(j-i)];
        }
        s << endl;
    }
}

// convert a symmetric matrix A to a full matrix B, such that mathematically A == B, where
//    A is the nr-by-nc input symmetric matrix stored in sa[]
//    B is the nr-by-nc output general matrix stored in sb[]
// if sb==NULL is given, then the required memory of size nr*nc will be allocated for the output matrix B;
//    otherwise it is assumed that sufficient memory has been allocated, with the address pointed to by sb
void symmetric_to_general(mkIndex nrc, const Real *sa, Real *&sb) {
    if (nrc == 0) {
        // a special case, empty input, do nothing
        return;
    }

    if (sb == NULL) {
        // allocate memory
        sb = new Real[nrc*nrc];
    }

    // now convert symmetric A to general B
    const Real *a0 = sa;
    Real *b0 = sb;
    for (mkIndex j=1; j<=nrc; j++) {
        // another pointer to B(j,j)
        Real *b2 = b0;

        // OCTAVE/MATLAB: B(j,j) = A(j,j)
        *b0++ = *a0++;

        for (mkIndex i=j+1; i<=nrc; i++) {
            // OCTAVE/MATLAB: B(i,j) = A(i,j)
            *b0++ = *a0;

            // OCTAVE/MATLAB: B(j,i) = A(i,j)
            b2 += nrc;
            *b2 = *a0++;
        }
        // for pointer to B(j+1,j+1)
        b0 += j;
    }
}

// compute a vector v formed by the k-th diagonal of a symmetric matrix A
// in one word, this routine computes v(i)=A(i+|k|,i) for i=1,...,nrc-|k|, where
//    A is the nrc-by-nrc input symmetric matrix stored in packed form in sa[], and
//    v is the output vector stored in sv[]
// the return value is the length of v
// if sv==NULL is given, then the required memory will be allocated for output vector v;
//    otherwise, it is assumed that sufficient memory has been allocated, with the address pointed to by sv
// this routine resembles v = diag(A,k) in OCTAVE/MATLAB
mkIndex symmatrix_diag(mkIndex nrc, const Real *sa, Real *&sv, mkSignedIndex k) {
    // absolute value of k
    mkIndex k0 = (k>=0) ? k : -k;

    if (nrc <= k0) {
        // a special case, the output vector is empty
        return 0;
    }

    // length of the output vector
    mkIndex len = nrc-k0;

    if (sv == NULL) {
        // allocate memory
        sv = new Real[len];
    }

    // for address of A(1+|k|,1)
    sa += k0;
    for (mkIndex i=0; i<len; i++) {
        // copy A(i+1+|k|,i+1) to sv[i]
        sv[i] = *sa;
        // for the next address
        sa += (nrc-i);
    }
    // else no element to copy, so do nothing

    // return the length of the vector
    return len;
}

// compute a vector v formed by the elements i,...,j of row/column k of an nrc-by-nrc symmetric matrix A
// 1. the input is an nrc-by-nrc symmetric matrix A stored in packed form in sa[]
// 2. the output is v=A(i:j,k) if indexFrom1==true, or v=A(i+1:j+1,k+1) if indexFrom1==false,
//    where v is a vector stored in sv[0],sv[incv],...,sv[|j-i|*incv]
// the indices i,j,k range over 1,...,nrc (if indexFrom1==true) , or over 0,...,nrc-1 (if indexFrom1==false)
// remarks:
// 1. if sv==NULL is given, then the required memory of |j-i|+1 will be allocated for the output vector v, and
//    incv will be set as 1 regardless of its input value;
//    otherwise it is assumed that all memory accesses via the pointer sv are legitimate
// 2. this routine allows i>j, in which case
//    the output is v=A(i:-1:j,k) if indexFrom1==true, or B=A(i+1:-1:j+1,k+1) if indexFrom1==false
void elements_of_symmatrix(mkIndex nrc, const Real *sa, mkIndex k, mkIndex i, mkIndex j, Real *&sv, mkSignedIndex incv, bool indexFrom1) {
    if (indexFrom1) {
        // shift the indices to start from 0
        i --;  j --;  k --;
    }

    // set flags according to whether i<=j or not
    int step;
    mkIndex ii, jj;
    if (i <= j) {
        ii = i;  jj = j;  step = 1;
    }
    else {
        ii = j;  jj = i;  step = -1;
    }

    // length of output
    mkIndex len = (jj-ii)+1u;

    if (sv == NULL) {
        // allocate memory
        sv = new Real[len];
        incv = 1;
    }

    // find the address to which A(ii+1,k+1) is copied
    Real *v0 = sv;
    if (step == -1) {
        v0 += (len-1)*incv;
        incv *= -1;
    }

    // find the address of A(k+1,ii+1) if ii<k, or A(ii+1,k+1) if ii>=k
    // if ii<k, also copy the elements A(k+1,ii+1:k)
    const Real *a0 = sa;
    if (ii < k) {
        a0 += ((ii*(2u*nrc-ii+1)/2u) + (k-ii));
        for (mkIndex z=ii; z<k; z++) {
            *v0 = *a0;
            a0 = a0 + nrc - z - 1;
            v0 += incv;
        }
        len -= (k-ii);
    }
    else {
        a0 += ((k*(2u*nrc-k+1)/2u) + (ii-k));
    }

    // if ii>=k, copy all A(ii:jj,k)
    // if ii<k, copy the rest elements, A(k:jj:,k)
    while (len--) {
        *v0 = *a0++;
        v0 += incv;
    }
}

// compute a general matrix B as the submatrix formed by elements in rows i1,...,i2 and columns j1,...,j2 of a symmetric matrix A
// 1. the input is an nrc-by-nrc symmetric matrix A stored in packed form in sa[]
// 2. the output is B=A(i1:i2,j1:j2) if indexFrom1==true, or B=A(i1+1:i2+1,j1+1:j2+1) if indexFrom1==false,
//    where B is a general matrix stored in sb[]
// the indices i1,i2,j1,j2 range over 1,...,nrc (if indexFrom1==true), or over 0,...,nrc-1 (if indexFrom1==false)
// remarks:
// 1. if sb==NULL is given, then the required memory of size (|i2-i1|+1)*(|j2-j1|+1) will be allocated for the output matrix B;
//    otherwise it is assumed that sufficient memory has been allocated, with the address pointed to by sb
// 2. this routine allows i1>i2 and/or j1>j2
//    if i1>i2 and j2>j2, then the output is B=A(i1:-1:i2,j1:-1:j2) if indexFrom1==true, or B=A(i1+1:-1:i2+1,j1+1:-1:j2+1) if indexFrom1==false
void submatrix_of_symmatrix(mkIndex nrc, const Real *sa, mkIndex i1, mkIndex i2, mkIndex j1, mkIndex j2, Real *&sb, bool indexFrom1) {
    // find the number of rows and the number of columns of the submatrix
    mkIndex nr = (i1<=i2) ? (i2-i1+1u) : (i1-i2+1u);
    mkIndex nc = (j1<=j2) ? (j2-j1+1u) : (j1-j2+1u);

    if (sb == NULL) {
        // allocate memory
        sb = new Real[nr*nc];
    }

    Real *b0 = sb;
    if (nr >= nc) {
        // copy column-by-column
        mkIndex j = j1;
        while (j != j2) {
            elements_of_symmatrix(nrc, sa, j, i1, i2, b0, 1u, indexFrom1);
            b0 += nr;
            if (j1 < j2)
                j ++;
            else
                j --;
        }
        // copy the last column
        elements_of_symmatrix(nrc, sa, j, i1, i2, b0, 1u, indexFrom1);
    }
    else {
        // copy row-by-row
        mkIndex i = i1;
        while (i != i2) {
            elements_of_symmatrix(nrc, sa, i, j1, j2, b0, nr, indexFrom1);
            b0 ++;
            if (i1 < i2)
                i ++;
            else
                i --;
        }
        // copy the last row
        elements_of_symmatrix(nrc, sa, i, j1, j2, b0, nr, indexFrom1);
    }
}

// compute a symmetric matrix B as the submatrix formed by elements in rows k1,...,k2 and columns k1,...,k2 of a symmetric matrix A
// 1. the input is an nrc-by-nrc symmetric matrix A stored in packed form in sa[]
// 2. the output is B=A(k1:k2,k1:k2) if indexFrom1==true, or B=A(k1+1:k2+1,k1+1:k2+1) if indexFrom1==false,
//    where B is a symmetric matrix stored in packed form in sb[]
// the indices k1,k2 range over 1,...,nrc (if indexFrom1==true), or over 0,...,nrc-1 (if indexFrom1==false)
// remarks:
// 1. if sb==NULL is given, then the required memory of size (|k2-k1|+1)*(|k2-k1|+2)/2 will be allocated for the output symmetric matrix B;
//    otherwise it is assumed that sufficient memory has been allocated, with the address pointed to by sb
// 2. this routine allows k1>k2, in which case
//    the output is B=A(k1:-1:k2,k1:-1:k2) if indexFrom1==true, or B=A(k1+1:-1:k2+1,k1+1:-1:k2+1) if indexFrom1==false
void submatrix_of_symmatrix(mkIndex nrc, const Real *sa, mkIndex k1, mkIndex k2, Real *&sb, bool indexFrom1) {
    // find the dimension of the submatrix
    mkIndex nrcb = (k1<=k2) ? (k2-k1+1u) : (k1-k2+1u);

    if (indexFrom1) {
        // shift the indices to start from 0
        k1 --;  k2 --;
    }

    if (sb == NULL) {
        // allocate memory
        sb = new Real[nrcb*(nrcb+1u)/2u];
    }

    if (k1 <= k2) {
        // for address of A(k1+1,k1+1)
        const Real *a0 = sa + k1*(2u*nrc-k1+1u)/2u;
        // for address of B(1,1)
        Real *b0 = sb;
        mkIndex stride = nrc-1-k2;
        for (mkIndex j=k1; j<=k2; j++) {
            // copy A(j+1:k2+1,j+1) to B(j-k1+1:nrcb,j-k1+1)
            mkIndex len = k2-j+1;
            for (mkIndex i=0; i<len; i++) {
                // copy A(j+1+i,j+1) to B(j-k1+1+i,j-k1+1)
                *b0++ = *a0++;
            }
            // for address of A(j+2,j+2) for the next iteration
            a0 += stride;
        }
    }
    else {  // k2 < k1
        // for address of A(k2+1,k2+1)
        const Real *a0 = sa + k2*(2u*nrc-k2+1u)/2u;
        mkIndex stride = nrc-1-k1;
        for (mkIndex j=nrcb; j>0; j--) {
            // for address of B(j,j)
            Real *b0 = sb + (j-1u)*(2*nrcb-j+2u)/2u;
            // copy A(k2+1+nrcb-j:k2+nrcb,k2+1+nrcb-j) to B(j,j:-1:1)
            for (mkIndex i=1; i<=j; i++) {
                // copy A(k2+nrcb-j+i,k2+1+nrcb-j) to B(j,j-i+1)
                *b0 = *a0++;
                // for address of B(j,j-i)
                b0 -= i;
            }
            // for address of A(k2+nrcb-j+2,k2+nrcb-j+2) for the next iteration
            a0 += stride;
        }
    }
}

// symmetrically permute rows and columns of a symmetric matrix A (not in-place)
// the result is a symmetric matrix B, such that
//    B(perm[i-1],  perm[j-1]  ) = A(i,j) for i,j=1,...,nrc (if indexFrom1==true), or
//    B(perm[i-1]+1,perm[j-1]+1) = A(i,j) for i,j=1,...,nrc (if indexFrom1==false), where
//    A is the  input symmetric matrix stored in packed form in sa[], and
//    B is the output symmetric matrix stored in packed form in sb[]
// the indices perm[0,...,nrc-1] range over 1,...,nrc (if indexFrom1==true), or over 0,...,nrc-1 (if indexFrom1==false)
// if sb==NULL is given, then the required memory of size nrc*(nrc+1)/2 will be allocated for the output symmetric matrix B;
//    otherwise it is assumed that sufficient memory has been allocated, with the address pointed to by sb
// this routine resembles B(perm,perm) = A in OCTAVE/MATLAB
void permute_symmatrix_rows_and_columns(mkIndex nrc, const mkIndex *perm, const Real *sa, Real *&sb, bool indexFrom1) {
    if (nrc == 0) {
        // a special case, empty input, do nothing
        return;
    }

    if (sb == NULL) {
        // allocate memory
        sb = new Real[nrc*(nrc+1u)/2u];
    }

    for (mkIndex j=0; j<nrc; j++) {
        // copy A(j+1:nrc,j+1) to B(perm[j+1:nrc],perm[j]) if indexFrom1==true, or
        // copy A(j+1:nrc,j+1) to B(perm[j+1:nrc]-1,perm[j]-1) if indexFrom1==false
        mkIndex j2 = perm[j];
        if (indexFrom1)
            j2 --;
        for (mkIndex i=j; i<nrc; i++) {
            mkIndex i2 = perm[i];
            if (indexFrom1)
                i2 --;
            // copy A(i+1,j+1) to B(i2+1,j2+1)
            mkIndex idx2 = (i2 >= j2) ?
                           ((i2-j2) + ((2u*nrc-j2+1u)*j2)/2u) :
                           ((j2-i2) + ((2u*nrc-i2+1u)*i2)/2u);
            sb[idx2] = *sa++;
        }
    }
}



////////////////////////////////////////////////////////////////////////////
//    sparse matrix routines
////////////////////////////////////////////////////////////////////////////

// print out sparse matrix A to ostream s, where A is of size nr-by-nc in CSC format stored in sa[], rowIdx[], colPtr[]
// if patternOnly==true, then only the sparsity pattern is to be printed, and sa[] will not be used
// if indexFrom1==true (or indexFrom1==false), the printed row/column indices start from 1 (or 0, respectively)
void print_spmatrix(std::ostream &s, mkIndex nr, mkIndex nc, const Real *sa, const mkIndex *rowIdx, const mkIndex *colPtr, bool indexFrom1, bool patternOnly) {
    std::ostream *sptr = &s;
    #ifdef USE_MEX
    if (&s == &std::cout || &s == &std::cerr) {
        // stdout and stderr are not supported by the mex interface
        // use mexostream (which inherits ostream) "mexcout" for printing messages to OCTAVE/MATLAB console
        // see matkitdef.h for more information
        sptr = &mexcout;
    }
    #endif

    *sptr << "compressed sparse columns (rows = " << nr << ", columns = " << nc << ", nnz = " << colPtr[nc]-colPtr[0] << ")" << endl;
    for (mkIndex j=0; j<nc; j++) {
        mkIndex z = colPtr[j+1]-colPtr[0];
        for (mkIndex ii=colPtr[j]-colPtr[0]; ii<z; ii++) {
            *sptr << "    (" << rowIdx[ii]+indexFrom1 << "," << j+indexFrom1 << ")";
            if (!patternOnly)
                *sptr << "     " << sa[ii];
            *sptr << endl;
        }
    }
}

// write a sparse matrix A to a matrix-market file fileName[],
//    where A is in CSC format stored in store[], rowIdx[], colPtr[]
// nrows, ncols are the number of rows, the number of columns, respectively
// if patternOnly==true, then only the sparsity pattern is to be written to the file fileName[], and
//    sa[] will not be used
// comment[] is a character string which will be printed in the output file fileName[] (if comment!=NULL)
// the return value:
//    0: a successful write
//    1: file open error
//    2: file close error
//    a negative integer: a number returned by fprintf() which signals an error
int matrix_market_spmatrix_write(mkIndex nrows, mkIndex ncols, const Real *sa, const mkIndex *rowIdx, const mkIndex *colPtr, const char fileName[], bool patternOnly, const char comment[]) {
    // open the file to write
    FILE *mmFile = fopen(fileName, "w");
    if (mmFile == NULL) {
        // file open error!
        return 1;
    }

    // print the header
    int info;
    if (patternOnly)
        info = fprintf(mmFile, "%%%%MatrixMarket matrix pattern real general\n");
    else
        info = fprintf(mmFile, "%%%%MatrixMarket matrix coordinate real general\n");
    if (info < 0) {
        // file write error
        fclose(mmFile);
        return info;
    }

    // print the comment
    info = fprintf(mmFile, "%%-------------------------------------------------------------------------------\n");
    if (info >= 0) {
        if (comment == NULL)
            info = fprintf(mmFile, "%% matrix-market sparse matrix file generated by matrix_market_spmatrix_write()\n");
        else
            info = fprintf(mmFile, "%% %s\n", comment);
    }
    if (info >= 0)
        info = fprintf(mmFile, "%%-------------------------------------------------------------------------------\n");
    if (info < 0) {
        // file write error
        fclose(mmFile);
        return info;
    }

    // print nrows (number of rows), ncols (number of columns), and nnz (number of nonzero elements)
    fprintf(mmFile, "%u %u %u\n", (unsigned)nrows, (unsigned)ncols, (unsigned)(colPtr[ncols]-colPtr[0]));

    // print the sparse matrix
    for (mkIndex j=1; j<=ncols; j++) {
        mkIndex z = colPtr[j]-colPtr[0];
        for (mkIndex ii=colPtr[j-1u]-colPtr[0]; ii<z; ii++) {
            if (patternOnly) {
                info = fprintf(mmFile, "%u %u\n", (unsigned)rowIdx[ii]+1u, j);
            }
            else {
                #ifdef USE_SINGLE
                info = fprintf(mmFile, "%u %u %.7g\n",   (unsigned)rowIdx[ii]+1u, (unsigned)j, sa[ii]);
                #else
                info = fprintf(mmFile, "%u %u %.15lg\n", (unsigned)rowIdx[ii]+1u, (unsigned)j, sa[ii]);
                #endif
            }
            if (info < 0) {
                // file write error
                fclose(mmFile);
                return info;
            }
        }
    }

    // now close the file
    if (fclose(mmFile) == EOF) {
        // close file error
        return 2;
    }

    // a successful write, return 0
    return 0;
}

// read a sparse matrix from a matrix-market file fileName[] and store the result in the coordinate format in sa[], ridx[], cidx[]
// the number of rows, the number of columns are recorded as nrows, ncols, respectively
// symmetry information will also be recorded, with symm='g' for general, symm='s' for symmetric, and symm='k' for skew-symmetric
//    for the latter two cases ('s' and 'k'), each pair of off-diagonal elements A(i,j) and A(j,i) is recorded once in the output
// the return value:
//    0: a successful read
//    1: the input file is a dense matrix
//    2: the input file is a complex (sparse) matrix
//    3: the input file contains only sparsity information
//   -1: file open error
//   -2: file is empty
//   -3: the first line (the header) of the file is empty
//   -4: invalid header
//   -5: invalid number of rows, number of columns, or number of nonzero elements
//   -6: the matrix data is incomplete or invalid
// if sa==NULL (ridx==NULL, or cidx==NULL) is given, then the required memory will be allocated;
//    otherwise it is assumed that sufficient memory has been allocated, with the address pointed to by sa (ridx, cidx, respectively)
// it also takes matrix files in the straight coordinate format (i.e., without the matrix-market header), in which case symm='g' will be set
int matrix_market_spmatrix_read(const char fileName[], mkIndex &nrows, mkIndex &ncols, mkIndex &nnz, Real *&sa, mkIndex *&ridx, mkIndex *&cidx, char &symm) {
    // open the file to read
    FILE *mmFile = fopen(fileName, "r");
    if (mmFile == NULL) {
        // file open error!
        return -1;
    }

    char string[256];  // assume that there are at most 255 characters in each line of the header
    if (fgets(string, 256, mmFile) == NULL) {
        // file is empty!
        fclose(mmFile);
        return -2;
    }

    if (strlen(string) == 0) {
        // the first line in the input file is empty!
        fclose(mmFile);
        return -3;
    }

    if (*string != '%') {
        // assume that the file has flat coordinates and values
        symm = 'g';
    }
    else {
        char *token = strtok(string, " \n");  // `%%MatrixMarket' should be read
        if (token == NULL || strcmp(token, "%%MatrixMarket")) {
            // invalid header
            fclose(mmFile);
            return -4;
        }

        token = strtok(NULL, " \n");  // "matrix" should be read
        if (token == NULL) {
            // invalid header
            fclose(mmFile);
            return -4;
        }
        size_t len = strlen(token);
        for (mkIndex i=0; i<len; i++)
            token[i] = tolower(token[i]);
        if (strcmp(token, "matrix")) {
            // not "matrix"!?  invalid header
            fclose(mmFile);
            return -4;
        }

        token = strtok(NULL, " \n");  // "coordinate" should be read, for a sparse matrix
                                      // the other valid matrix-market tag is "array" for a general matrix
        if (token == NULL) {
            // invalid header
            fclose(mmFile);
            return -4;
        }
        len = strlen(token);
        for (mkIndex i=0; i<len; i++)
            token[i] = tolower(token[i]);
        if (!strcmp(token, "array")) {
            // input matrix-market file stores a dense matrix; skip the rest
            fclose(mmFile);
            return 1;
        }
        if (strcmp(token, "coordinate")) {
            // not "coordinate"!?  invalid header
            fclose(mmFile);
            return -4;
        }

        token = strtok(NULL, " \n");  // "real" should be read, for a real-valued matrix, "integer" is also accepted
                                      // the other valid matrix-market tags are "complex" and "pattern"
        if (token == NULL) {
            // invalid header
            fclose(mmFile);
            return -4;
        }
        len = strlen(token);
        for (mkIndex i=0; i<len; i++)
            token[i] = tolower(token[i]);
        if (!strcmp(token, "complex")) {
            // the input matrix-market file stores a complex (sparse) matrix
            fclose(mmFile);
            return 2;
        }
        if (!strcmp(token, "pattern")) {
            // the input matrix-market file contains only the sparsity information
            fclose(mmFile);
            return 3;
        }
        if (strcmp(token, "real") && strcmp(token, "integer")) {
            // neither "real" nor "integer"!?  invalid header
            fclose(mmFile);
            return -4;
        }

        token = strtok(NULL, " \n");  // "general", "symmetric", "skew-symmetric" should be read, "hermitian" is also possible
                                      // "skew-hermitian" seems not a valid tag in the matrix-market format, but is allowed here
                                      // however, "hermitian" and "skew-hermitian" should not appear here, since we are reading a real-valued matrix
        if (token == NULL) {
            // invalid header
            fclose(mmFile);
            return -4;
        }
        len = strlen(token);
        for (mkIndex i=0; i<len; i++)
            token[i] = tolower(token[i]);
        if (!strcmp(token, "general")) {
            // input matrix-market file stores a general matrix
            symm = 'g';
        }
        else if (!strcmp(token, "symmetric") || !strcmp(token, "hermitian")) {
            // input matrix-market file stores a symmetric matrix
            symm = 's';
        }
        else if (!strcmp(token, "skew-symmetric") || !strcmp(token, "skew-hermitian")) {
            // input matrix-market file stores a symmetric matrix
            symm = 'k';
        }
        else {
            // not a recognized tag;  invalid header
            fclose(mmFile);
            return -4;
        }

        // get rid of the rest comment lines
        do {
            if (fgets(string, 256, mmFile) == NULL) {
                // no information about nrows, ncols, nnz!?
                fclose(mmFile);
                return -5;
            }
        } while (strlen(string)==0 || *string=='\n' || *string=='%');
    }

    // now parse the line which contains nrows, ncols, nnz
    // read nrows
    char *token = strtok(string, " \n");
    if (token == NULL) {
        // no information about the number of rows!?
        fclose(mmFile);
        return -5;
    }
    int i = atoi(token);
    if (i <= 0) {
        // number of rows is not positive!?
        fclose(mmFile);
        return -5;
    }
    nrows = (mkIndex)i;

    // read ncols
    token = strtok(NULL, " \n");
    if (token == NULL) {
        // no information about the number of columns!?
        fclose(mmFile);
        return -5;
    }
    i = atoi(token);
    if (i <= 0) {
        // number of columns is not positive!?
        fclose(mmFile);
        return -5;
    }
    ncols = (mkIndex)i;

    // read number of nonzero elements (nnz)
    token = strtok(NULL, " \n");
    if (token == NULL) {
        // no information about nnz!?
        fclose(mmFile);
        return -5;
    }
    mkSignedIndex nnz0 = atoi(token);
    if (nnz0 < 0) {
        // nnz is negative!?
        fclose(mmFile);
        return -5;
    }
    nnz = nnz0;

    // allocate memory if necessary
    if (nnz) {
        if (ridx == NULL)
            ridx = new mkIndex[nnz];
        if (cidx == NULL)
            cidx = new mkIndex[nnz];
        if (sa == NULL)
            sa = new Real[nnz];
    }

    // read the matrix
    mkIndex *r = ridx;
    mkIndex *c = cidx;
    Real *s = sa;
    while (nnz0--) {
        #ifdef USE_SINGLE
        if (fscanf(mmFile, "%u %u %g\n",  r++, c++, s++) != 3) {
        #else
        if (fscanf(mmFile, "%u %u %lg\n", r++, c++, s++) != 3) {
        #endif
            // incomplete data!?
            fclose(mmFile);
            return -6;
        }
    }

    // a successful read, close the input file and return 0
    fclose(mmFile);
    return 0;
}

// convert a sparse matrix in coordinate format to the compressed sparse column (CSC) format
// key parameters:
// 1. nc and nnz (as input) are number of columns and number of nonzero elements
// 2. the input coordinate format is specified by sa[], ridx[], cidx[]
// 3. the output CSC format is specified by store[], rowIdx[], colPtr[]
// return:
//    0 on success, or -1 if symm is not valid
// other parameters:
// 1. symm gives the symmetry information: 'g' for general (default), 's' for symmetric, and 'k' for skew-symmetric
//    if A is symmetric or skew-symmetric, then each pair of nonzero off-diagonal elements A(i,j) and A(j,i) is recorded once in the coordinate format (as input);
//    however the elements will be duplicated in the output CSC format (as output) for being symmetric or skew-symmetric
// 2. if indexFrom1==true (default),  the row indices ridx[] and column indices cidx[] range over 1,...,nr   and 1,...,nc,   respectively;
//    if indexFrom1==false,           the row indices ridx[] and column indices cidx[] range over 0,...,nr-1 and 0,...,nc-1, respectively
//    the output row indices rowIdx[] always range over 0,...,nr-1
// 3. if patternOnly==true (default is false), then only the sparsity pattern is considered, and sa[] and store[] will not be used
// remarks:
// 1. the number of rows of A, nr, is not required in the computation, so it is not passed
// 2. if store==NULL (rowIdx==NULL, or colPtr==NULL) is given, then the required memory will be allocated;
//    otherwise it is assumed that sufficient memory has been allocated, with the address pointed to by store (rowIdx, or colPtr, respectively)
// 3. this routine can also convert a sparse matrix in coordinate format to the compressed sparse row (CSR) format,
//    by swapping ridx[] and cidx[] and passing nc as the number of rows
int coordinates_to_csc(mkIndex nc, mkIndex nnz0, const Real *sa, const mkIndex *ridx, const mkIndex *cidx, Real *&store, mkIndex *&rowIdx, mkIndex *&colPtr, char symm, bool indexFrom1, bool patternOnly) {
    // find the `true' number of nonzero elements, assuming all recorded elements are non-zero, and
    // if A is symmetric or skew-symmetric, then each pair of nonzero off-diagonal elements A(i,j) and A(j,i) is recorded once in the coordinate format
    mkIndex nnz;
    if (symm == 'g' || symm == 'G') {
        // general matrix
        nnz = nnz0;
    }
    else if (symm == 's' || symm == 'S') {
        // symmetric matrix
        nnz = 2*nnz0;
        for (mkIndex z=0; z<nnz0; z++) {
            if (ridx[z] == cidx[z])
                nnz --;
        }
    }
    else if (symm == 'k' || symm == 'K') {
        // skew-symmetric matrix
        nnz = 0;
        for (mkIndex z=0; z<nnz0; z++) {
            if (ridx[z] != cidx[z])
                nnz ++;
        }
        nnz *= 2;
    }
    else {
        // Error: symm must be 'g' (or 'G') for general, 's' (or 'S') for symmetric, or 'k' (or 'K') for skew-symmetric!
        return -1;
    }

    // allocate memory for output, if not provided
    if (!patternOnly && store == NULL && nnz)
        store = new Real[nnz];
    if (rowIdx == NULL && nnz)
        rowIdx = new mkIndex[nnz];
    if (colPtr == NULL)
        colPtr = new mkIndex[nc+1u];

    // zero initialization of colPtr[0,...,nc]
    for (mkIndex j=0; j<=nc; j++)
        colPtr[j] = 0;

    // find the number of nonzero elements for each column
    mkIndex *cptr;
    if (indexFrom1)
        cptr = colPtr;
    else
        cptr = colPtr + 1;
    if (symm == 'g' || symm == 'G') {
        // general matrix
        for (mkIndex z=0; z<nnz0; z++)
            cptr[cidx[z]] ++;
    }
    else if (symm == 's' || symm == 'S') {
        // symmetric matrix
        for (mkIndex z=0; z<nnz0; z++) {
            cptr[cidx[z]] ++;
            if (cidx[z] != ridx[z])
                cptr[ridx[z]] ++;
        }
    }
    else {
        // skew-symmetric matrix (symm == 'k' || symm == 'K')
        for (mkIndex z=0; z<nnz0; z++) {
            if (cidx[z] != ridx[z]) {
                cptr[cidx[z]] ++;
                cptr[ridx[z]] ++;
            }
        }
    }

    // get the pointers to the positions
    for (mkIndex j=1; j<nc; j++)
        colPtr[j+1u] += colPtr[j];

    // fill the element values store[0,...,nnz-1] (if !patternOnly) and row indices rowIdx[0,...,nnz-1] of the matrix
    while (nnz0--) {
        Real a = *sa++;
        mkIndex r = *ridx++;
        mkIndex c = *cidx++;
        if (indexFrom1) {
            r--;
            c--;
        }
        if ((symm!='k' && symm!='K') || r!=c) {
            // not skew-symmetric or off-diagonal
            mkIndex loc = colPtr[c]++;
            if (!patternOnly)
                store[loc] = a;   // copy the value
            rowIdx[loc] = r;      // copy the row index
        }
        if ((symm=='s' || symm=='S') && r!=c) {
            // symmetric and off-diagonal
            mkIndex loc = colPtr[r]++;
            if (!patternOnly)
                store[loc] = a;   // copy the value
            rowIdx[loc] = c;      // copy the row index
        }
        if ((symm=='k' || symm=='K') && r!=c) {
            // skew-symmetric and off-diagonal
            mkIndex loc = colPtr[r]++;
            if (!patternOnly)
                store[loc] = -a;  // copy the value
            rowIdx[loc] = c;      // copy the row index
        }
    }

    // the column pointers colPtr[0,...,nc] have been shifted; re-shift the values back
    for (mkIndex j=nc; j>0; j--)
        colPtr[j] = colPtr[j-1u];
    colPtr[0] = 0;

    // successful conversion
    return 0;
}

// sort elements in each column of a sparse matrix A in CSC format (in-place implementation), where
//    A is the nr-by-nc input sparse matrix in CSC format stored in store[], rowIdx[], colPtr[],
//    but the elements in each column can be in any order (i.e. unsorted), and
// the output is the same matrix in CSC format, using the input storage store[], rowIdx[], colPtr[],
//    with the elements in each column sorted with respect to the row indices
// return 0 on success, or -1 if the computed permutation is not valid
// if patternOnly==true, then only the sparsity pattern is considered, and store[] will not be used
// if work==NULL is given, then work space of size nnz+max(nnz,nr+1) will be allocated for sorting, and the space will be freed when sorting is done;
//    otherwise, it is assumed that sufficient memory has been allocated with the address pointed to by work, and it will not be freed in this routine
// this routine can be used to sort elements in each row of a sparse matrix in CSR format
int sort_csc_elements_in_each_column(mkIndex nr, mkIndex nc, Real *store, mkIndex *rowIdx, mkIndex *colPtr,
                                     bool patternOnly, mkIndex *work) {
    // allocate work space, if it is not given
    mkIndex nnz = colPtr[nc]-colPtr[0];
    bool wgiven = true;  // whether the work space is given or not
    if (work == NULL) {
        wgiven = false;
        mkIndex wsize = nnz + ((nnz>nr) ? nnz : (nr+1u));
        work = new mkIndex[wsize];
    }
    mkIndex *work2 = work + nnz;

    // zero initialization of work2[0,...,nr]
    mkIndex i = nr+1;
    mkIndex *wptr = work2;
    while (i--)
        *wptr++ = 0;

    // count number of nonzero elements in each row
    mkIndex *ridx = rowIdx;
    wptr = work2+1;
    mkIndex z = nnz;
    while (z--)
        wptr[*ridx++]++;

    // compute the row pointers as those in CSR format
    for (i=1; i<nr; i++)
        work2[i+1u] += work2[i];

    // get the indices work[0,...,nnz-1] in order as the order of elements in CSR format
    for (z=0; z<nnz; z++)
        work[work2[rowIdx[z]]++] = z;

    // find column indices work2[0,...,nnz-1], work2[i] is the column index of i-th element in CSC format
    for (mkIndex j=0; j<nc; j++) {
        for (z=colPtr[j]-colPtr[0]; z<colPtr[j+1u]-colPtr[0]; z++)
            work2[z] = j;
    }

    // compute the permutation work2[0,...,nnz-1] for elements in each column sorted with respect to row indices
    for (z=0; z<nnz; z++) {
        mkIndex idx = work[z];     // consider the idx-th element in CSC format
        mkIndex col = work2[idx];  // the column containing this element
        work2[idx] = colPtr[col]++;      // work2[idx] is the index/order of this element in CSC format
                                         // with elements in each column sorted with respect to the row indices
    }

    // the column pointers colPtr[] have been shifted; re-shift the values back
    for (mkIndex j=nc; j>0; j--)
        colPtr[j] = colPtr[j-1u];
    colPtr[0] = 0;

    // now perform the permutation
    if (in_place_permute_elements(nnz, work2, rowIdx, false) == -1) {
        // the computed permutation is not valid
        return -1;
    }
    if (!patternOnly && in_place_permute_elements(nnz, work2, store, false) == -1) {
        // the computed permutation is not valid
        return -1;
    }

    // free the work space, if it is not from the input
    if (!wgiven) {
        delete [] work;
    }

    // successful sort sort
    return 0;
}

// convert a sparse matrix A to a full matrix B, such that mathematically A == B, where
//    A is the nr-by-nc input sparse matrix stored in sa[], rowIdx[], colPtr[]
//    B is the nr-by-nc output general matrix stored in sb[]
// if sb==NULL is given, then the required memory of size nr*nc will be allocated for the output matrix B;
//    otherwise it is assumed that sufficient memory has been allocated, with the address pointed to by sb
// this routine resembles B = full(A) in OCTAVE/MATLAB
void sparse_to_full(mkIndex nr, mkIndex nc, const Real *sa, const mkIndex *rowIdx, const mkIndex *colPtr, Real *&sb) {
    if (nr == 0 || nc == 0) {
        // a special case, empty input, do nothing
        return;
    }

    // size of the storage
    mkIndex sz = nr*nc;

    if (sb == NULL) {
        // allocate memory
        sb = new Real[sz];
    }

    // zero initialization
    Real *b0 = sb;
    while (sz--)
        *b0++ = 0.0;

    // now convert sparse A to full B
    b0 = sb;
    const mkIndex *ridx = rowIdx;
    const Real *a0 = sa;
    for (mkIndex j=1; j<=nc; j++) {
        // OCTAVE/MATLAB: B(:,j) = A(:,j)
        mkIndex z = colPtr[j]-colPtr[j-1u];
        while (z--) {
            // OCTAVE/MATLAB: B(*ridx,j) = A(*ridx,j)
            b0[*ridx++] = *a0++;
        }
        // for pointer to B(:,j+1)
        b0 += nr;
    }
}

// compute a vector v formed by the k-th diagonal of a sparse matrix A
// to be precise, this routine computes v(i)=A(i,i+k) if k>=0, or v(i)=A(i-k,i) if k<0, where
//    A is the nr-by-nc input sparse matrix in CSC format stored in store[], rowIdx[], colPtr[]
//    v is the output vector v stored in sv[] containing the k-th diagonal elements (in order) of A
// the return value is the length of v
// if sv==NULL is given, then the required memory will be allocated for the output vector v;
//    otherwise, it is assumed that sufficient memory has been allocated, with the address pointed to by sv
// this routine can also be used to find the k-th diagonal of a sparse matrix in CSR format, by changing the sign of k
// this routine resembles v = diag(A,k) in OCTAVE/MATLAB
mkIndex spmatrix_diag(mkIndex nr, mkIndex nc, const Real *sa, const mkIndex *rowIdx, const mkIndex *colPtr, Real *&sv, mkSignedIndex k) {
    // compute len: the length of output vector
    //          k0: (k>=0) ? k : 0
    //          k1: (k>=0) ? 0 : -k
    mkIndex len, k0, k1;
    if (k >= 0) {
        len = (nc>(mkIndex)k) ? (nc-(mkIndex)k) : 0;
        if (len > nr)
            len = nr;
        k0 = k;
        k1 = 0;
    }
    else {
        len = (nr>(mkIndex)(-k)) ? (nr-(mkIndex)(-k)) : 0;
        if (len > nc)
            len = nc;
        k0 = 0;
        k1 = -k;
    }

    if (sv == NULL && len) {
        // allocate memory
        sv = new Real[len];
    }

    // copy the elements
    for (mkIndex j=0; j<len; j++) {
        sv[j] = 0.0;
        mkIndex z = colPtr[j+k0+1u]-colPtr[0];
        for (mkIndex ii=colPtr[j+k0]-colPtr[0]; ii<z; ii++) {
            if (rowIdx[ii] == j+k1)
                sv[j] = sa[ii];
            if (rowIdx[ii] > j+k1)
                break;  // break the inner loop since the (j+1)st diagonal element is zero, assuming that the row indices are sorted in each column
        }
    }

    // return the length of output vector
    return len;
}

// transpose a sparse matrix A and store the result as a sparse matrix B (not in-place)
// to be precise, mathematically B(j,i)=A(i,j) for i=1,...,nr and j=1,...,nc, where
//    A is the nr-by-nc input sparse matrix in CSC format stored in store[], rowIdx[], colPtr[]
//    B, the transpose of A, is the nc-by-nr output sparse matrix in CSC format stored in store2[], rowIdx2[], colPtr2[]
// if patternOnly==true, then only the sparsity pattern is considered, and store[] and store2[] will not be used
// if store2==NULL (rowIdx2==NULL, or colPtr2==NULL) is given, then the required memory will be allocated;
//    otherwise it is assumed that sufficient memory has been allocated, with the address pointed to by store2 (rowIdx2, or colPtr2, respectively)
// the output sparse matrix B has elements in each column sorted with respect to row indices,
//    no matter whether the input matrix has this property or not
// this routine resembles B = A' or equivalently B = transpose(A) in OCTAVE/MATLAB
// this routine can also be used for transposing a sparse matrix in CSR format
void spmatrix_transpose(mkIndex nr, mkIndex nc, const Real *store, const mkIndex *rowIdx, const mkIndex *colPtr,
                        Real *&store2, mkIndex *&rowIdx2, mkIndex *&colPtr2, bool patternOnly) {
    // allocate memory for the matrix transpose
    mkIndex nnz = colPtr[nc]-colPtr[0];
    if (!patternOnly && store2 == NULL && nnz)
        store2 = new Real[nnz];
    if (rowIdx2 == NULL && nnz)
        rowIdx2 = new mkIndex[nnz];
    if (colPtr2 == NULL)
        colPtr2 = new mkIndex[nr+1u];
    // the input matrix has nr rows, so its transpose, the output matrix has nr columns

    // zero initialization of colPtr2
    mkIndex j = nr+1u;
    mkIndex *cptr = colPtr2;
    while (j--)
        *cptr++ = 0;

    // count number of nonzero elements in each row of the input matrix
    const mkIndex *ridx = rowIdx;
    cptr = colPtr2+1;
    while (nnz--)
        cptr[*ridx++] ++;

    // compute column pointers of the output matrix
    for (j=1; j<nr; j++)
        colPtr2[j+1u] += colPtr2[j];

    // do the transpose
    for (mkIndex i=0; i<nc; i++) {
        for (mkIndex jj=colPtr[i]-colPtr[0]; jj<colPtr[i+1u]-colPtr[0]; jj++) {
            j = rowIdx[jj];
            mkIndex loc = colPtr2[j]++;
            if (!patternOnly)
                store2[loc] = store[jj];  // copy the value
            rowIdx2[loc] = i;  // copy the column index of the input matrix to the row index of the output matrix
        }
    }

    // the column pointers colPtr2[] of the output matrix have been shifted; re-shift the values back
    for (j=nr; j>0; j--)
        colPtr2[j] = colPtr2[j-1u];
    colPtr2[0] = 0;
}

// compute a sparse matrix B as the submatrix formed by elements in rows i1,...,i2 and columns j1,...,j2 of a sparse matrix A
// 1. the input is an nr-by-nc sparse matrix A stored in CSC format in store[], rowIdx[], colPtr[]
// 2. the output is B=A(i1:i2,j1:j2) if indexFrom1==true, or B=A(i1+1:i2+1,j1+1:j2+1) if indexFrom1==false,
//    where B is a sparse matrix in CSC format stored in store2[], rowIdx2[], colptr2[]
// other parameters:
// 1. the row indices i1,i2 and column indices j1,j2 range over 1,...,nr and 1,...,nc (if indexFrom1==true), or
//    over 0,...,nr-1 and 0,...,nc-1 (if indexFrom1==false), respectively
// 2. if patternOnly==true, then only sparsity pattern is considered, and store[] and store2[] will not be used
// remarks:
// 1. the number of rows, nr, is not required in the computation, so it is not passed
// 2. if store2==NULL (rowIdx2==NULL, or colPtr2==NULL) is given, then the required memory will be allocated;
//    otherwise it is assumed that sufficient memory has been allocated, with the address pointed to by store2 (rowIdx2, or colPtr2, respectively)
// 3. this routine allows i1>i2 and/or j1>j2
//    if i1>i2 and j2>j2, then the output is B=A(i1:-1:i2,j1:-1:j2) if indexFrom1==true, or B=A(i1+1:-1:i2+1,j1+1:-1:j2+1) if indexFrom1==false
// 4. this routine does not require the input sparse matrix with elements in each column sorted by the row indices
//    however the output is guaranteed to have this property, only if the input has this property
// 5. this routine can also be used to compute a submatrix of a sparse matrix in CSR format
void submatrix_of_spmatrix(mkIndex nc, const Real *store, const mkIndex *rowIdx, const mkIndex *colPtr,
                           mkIndex i1, mkIndex i2, mkIndex j1, mkIndex j2,
                           Real *&store2, mkIndex *&rowIdx2, mkIndex *&colPtr2, bool indexFrom1, bool patternOnly) {
    if (indexFrom1) {
        // working with indices starting from 0
        i1--;  i2--;  j1--;  j2--;
    }
    // number of columns of the submatrix
    mkIndex nc2 = (j1<=j2) ? (j2-j1+1) : (j1-j2+1);

    if (colPtr2 == NULL) {
        // allocate memory
        colPtr2 = new mkIndex[nc2+1u];
    }

    // compute the column pointers colPtr2[0,...,nc2] of the submatrix
    // zero initialization of colPtr2[]
    mkIndex j = nc2 + 1u;
    mkIndex *cptr2 = colPtr2;
    while (j--)
        *cptr2++ = 0;
    // compute number of nonzero elements in each column of the submatrix
    cptr2 = colPtr2 + 1;
    for (j=0; j<nc2; j++) {
        // consider B(:,j+1), the (j+1)st column of the submatrix
        mkIndex j0;
        if (j1 <= j2)
            j0 = j1 + j;
        else  // j1 > j2
            j0 = j1 - j;
        // A(:,j0+1) is the corresponding column of the original matrix
        mkIndex ii1 = (i1<=i2) ? i1 : i2;  // min(i1,i2)
        mkIndex ii2 = (i1<=i2) ? i2 : i1;  // max(i1,i2)
        mkIndex z = colPtr[j0+1] - colPtr[0];
        for (mkIndex k=colPtr[j0]-colPtr[0]; k<z; k++) {
            if (rowIdx[k]<=ii2 && rowIdx[k]>=ii1)
                cptr2[j] ++;
        }
    }
    // sum up numbers of nonzero elements to get column pointers
    for (j=1; j<nc2; j++)
        colPtr2[j+1u] += colPtr2[j];

    // number of nonzero elements of the submatrix
    mkIndex nnz2 = colPtr2[nc2];  // note that colPtr2[0] == 0 in all cases

    if (rowIdx2 == NULL && nnz2) {
        // allocate memory
        rowIdx2 = new mkIndex[nnz2];
    }
    if (!patternOnly && store2 == NULL && nnz2) {
        // allocate memory
        store2 = new Real[nnz2];
    }

    // set element values (if patternOnly==false) and row indices for the submatrix
    mkIndex *ridx2 = rowIdx2;
    Real *str2 = store2;
    for (mkIndex j=0; j<nc2; j++) {
        // set B(:,j+1), the (j+1)st column of the submatrix
        mkIndex j0;
        if (j1 <= j2)
            j0 = j1+j;
        else
            j0 = j1-j;
        // the corresponding column of the original matrix is the (j0+1)st column, A(:,j0+1)
        mkIndex z = colPtr[j0+1u] - colPtr[0];
        mkIndex k = colPtr[j0] - colPtr[0];
        if (i1 <= i2) {
            // compute B(1:i2-i1+1,j+1) = A(i1+1:i2+1,j0+1)
            while (k < z) {
                if (rowIdx[k]<=i2 && rowIdx[k]>=i1) {
                    *ridx2++ = rowIdx[k]-i1;
                    if (!patternOnly)
                        *str2++ = store[k];
                }
                k ++;
            }
        }
        else {
            // compute B(1:i1-i2+1,j+1) = A(i2+1:-1:i1+1,j0+1)
            while (z > k) {
                z --;
                if (rowIdx[z]<=i1 && rowIdx[z]>=i2) {
                    *ridx2++ = i1-rowIdx[z];
                    if (!patternOnly)
                        *str2++ = store[z];
                }
            }
        }
    }
}

// permute columns of a sparse matrix A in CSC format (not in-place)
// the result is a sparse matrix B, such that
//    B(i,perm[j-1])   = A(i,j) for i=1,...,nr and j=1,...,nc (if indexFrom1==true), or
//    B(i,perm[j-1]+1) = A(i,j) for i=1,...,nr and j=1,...,nc (if indexFrom1==false), where
//    A is the  input sparse matrix in CSC format stored in store[],  rowIdx[],  colPtr[], and
//    B is the output sparse matrix in CSC format stored in store2[], rowIdx2[], colPtr2[]
// other parameters:
// 1. the indices perm[0,...,nc-1] range over 1,...,nc (if indexFrom1==true), or over 0,...,nc-1 (if indexFrom1==false)
// 2. if patternOnly==true, then only the sparsity pattern is considered, and store[] and store2[] will not be used
// remarks:
// 1. the number of rows, nr, is not required in the computation, so it is not passed
// 2. if store2==NULL (rowIdx2==NULL, or colPtr2==NULL) is given, then the required memory will be allocated;
//    otherwise it is assumed that sufficient memory has been allocated, with the address pointed to by store2 (rowIdx2, or colPtr2, respectively)
// 3. this routine resembles B(:,perm) = A in OCTAVE/MATLAB
// 4. this routine can be used to permute rows of a sparse matrix in CSR format
void permute_csc_columns(mkIndex nc, const mkIndex *perm, const Real *store, const mkIndex *rowIdx, const mkIndex *colPtr,
                         Real *&store2, mkIndex *&rowIdx2, mkIndex *&colPtr2, bool indexFrom1, bool patternOnly) {
    // the column indices perm[0,...,nc-1] start from either 0 (indexFrom1==false) or 1 (indexFrom1==true)
    // that is, the values of perm[0,...,nc-1] are either 0,1,2,...,nc-1 (indexFrom1==false) or 1,2,3,...,nc (indexFrom1==true)

    // get nnz (number of nonzero elements)
    mkIndex nnz = colPtr[nc]-colPtr[0];

    // allocate memory for the output sparse matrix
    if (!patternOnly && store2 == NULL && nnz)
        store2 = new Real[nnz];
    if (rowIdx2 == NULL && nnz)
        rowIdx2 = new mkIndex[nnz];
    if (colPtr2 == NULL)
        colPtr2 = new mkIndex[nc+1u];

    // find the length of each column after permutation
    colPtr2[0] = 0;
    mkIndex *cptr2;
    if (indexFrom1)
        cptr2 = colPtr2;
    else
        cptr2 = colPtr2 + 1;
    for (mkIndex j=0; j<nc; j++)
        cptr2[perm[j]] = colPtr[j+1u]-colPtr[j];

    // get the right column pointers
    for (mkIndex j=0; j<nc; j++)
        colPtr2[j+1u] += colPtr2[j];

    // copy the permuted columns to the right locations
    cptr2 --;
    for (mkIndex j=0; j<nc; j++) {
        mkIndex len = colPtr[j+1u]-colPtr[j];  // column length
        mkIndex stride = colPtr[j]-colPtr[0];
        const Real *str = NULL;
        if (!patternOnly)
            str = store + stride;
        const mkIndex *ri = rowIdx + stride;
        Real *str2 = NULL;
        if (!patternOnly)
            str2 = store2 + cptr2[perm[j]];
        mkIndex *ri2 = rowIdx2 + cptr2[perm[j]];
        if (!patternOnly) {
            while (len--) {
                // copy elements
                *str2++ = *str++;
                // copy row indices
                *ri2++ = *ri++;
            }
        }
        else {
            while (len--)
                // copy row indices
                *ri2++ = *ri++;
        }
    }
}

// permute rows of a sparse matrix A in CSC format (not in-place)
// the result is a sparse matrix B, such that
//    B(perm[i-1],j)   = A(i,j) for i=1,...,nr and j=1,...,nc (if indexFrom1==true), or
//    B(perm[i-1]+1,j) = A(i,j) for i=1,...,nr and j=1,...,nc (if indexFrom1==false), where
//    A is the  input sparse matrix in CSC format stored in store[],  rowIdx[],  colPtr[], and
//    B is the output sparse matrix in CSC format stored in store2[], rowIdx2[], colPtr2[]
// other parameters:
// 1. the indices perm[0,...,nc-1] range over 1,...,nc (if indexFrom1==true), or over 0,...,nc-1 (if indexFrom1==false)
// 2. if patternOnly==true, then only the sparsity pattern is considered, and store[] and store2[] will not be used
// remarks:
// 1. the number of rows, nr, is not required in the computation, so it is not passed
// 2. if store2==NULL (rowIdx2==NULL, or colPtr2==NULL) is given, then the required memory will be allocated;
//    otherwise it is assumed that sufficient memory has been allocated, with the address pointed to by store2 (rowIdx2, or colPtr2, respectively)
// 3. the elements in each column of B may not be sorted with respect to the row indices, even if A has this property
//    if it is desired to have sorted elements in each column, then use sort_csc_elements_in_each_column() after permute_csc_rows()
// 4. this routine resembles B(perm,:) = A in OCTAVE/MATLAB
// 5. this routine can be used to permute columns of a sparse matrix in CSR format
void permute_csc_rows(mkIndex nc, const mkIndex *perm,
                      const Real *store, const mkIndex *rowIdx, const mkIndex *colPtr,
                      Real *&store2, mkIndex *&rowIdx2, mkIndex *&colPtr2, bool indexFrom1, bool patternOnly) {
    // the row indices perm[0,...,nr-1] start from either 0 (indexFrom1==false) or 1 (indexFrom1==true)
    // that is, the values of perm[0,...,nr-1] are either 0,1,2,...,nr-1 (indexFrom1==false) or 1,2,3,...,nr (indexFrom1==true)

    // get nnz (number of nonzero elements)
    mkIndex nnz = colPtr[nc]-colPtr[0];

    // allocate memory for the output sparse matrix
    if (!patternOnly && store2 == NULL && nnz)
        store2 = new Real[nnz];
    if (rowIdx2 == NULL && nnz)
        rowIdx2 = new mkIndex[nnz];
    if (colPtr2 == NULL)
        colPtr2 = new mkIndex[nc+1u];

    // permute the rows of the matrix A, and return the result
    if (!patternOnly)
        memory_xcopy(nnz, store, store2);
    memcpy(colPtr2, colPtr, (nc+1)*sizeof(mkIndex));
    const mkIndex *ri = rowIdx;
    mkIndex *ri2 = rowIdx2;
    if (indexFrom1) {
        while (nnz--)
            *ri2++ = perm[*ri++]-1u;
    }
    else {
        while (nnz--)
            *ri2++ = perm[*ri++];
    }
}

// permute rows and columns of a sparse matrix A in CSC format (not in-place)
// the result is a sparse matrix B, such that
//    B(rperm[i-1],cperm[j-1])   = A(i,j) for i=1,...,nr and j=1,...,nc (if indexFrom1==true), or
//    B(rperm[i-1],cperm[j-1]+1) = A(i,j) for i=1,...,nr and j=1,...,nc (if indexFrom1==false), where
//    A is the  input sparse matrix in CSC format stored in store[],  rowIdx[],  colPtr[], and
//    B is the output sparse matrix in CSC format stored in store2[], rowIdx2[], colPtr2[]
// other parameters:
// 1. the indices rperm[0,...,nr-1] range over 1,...,nr (if indexFrom1==true), or over 0,...,nr-1 (if indexFrom1==false)
// 2. the indices cperm[0,...,nc-1] range over 1,...,nc (if indexFrom1==true), or over 0,...,nc-1 (if indexFrom1==false)
// 3. if patternOnly==true, then only the sparsity pattern is considered, and store[] and store2[] will not be used
// remarks:
// 1. the number of rows, nr, is not required in the computation, so it is not passed
// 2. if store2==NULL (rowIdx2==NULL, or colPtr2==NULL) is given, then the required memory will be allocated;
//    otherwise it is assumed that sufficient memory has been allocated, with the address pointed to by store2 (rowIdx2, or colPtr2, respectively)
// 3. the elements in each column of B may not be sorted with respect to the row indices, even if A has this property
//    if it is desired to have sorted elements in each column, then use sort_csc_elements_in_each_column() after permute_csc_rows()
// 4. this routine resembles B(rperm,cperm) = A in OCTAVE/MATLAB
// 5. this routine can be used to permute columns and rows of a sparse matrix in CSR format
void permute_csc_rows_and_columns(mkIndex nc, const mkIndex *rperm, const mkIndex *cperm,
                                  const Real *store, const mkIndex *rowIdx, const mkIndex *colPtr,
                                  Real *&store2, mkIndex *&rowIdx2, mkIndex *&colPtr2, bool indexFrom1, bool patternOnly) {
    // the column indices cperm[0,...,nc-1] start from either 0 (indexFrom1==false) or 1 (indexFrom1==true)
    // that is, the values of cperm[0,...,nc-1] are either 0,1,2,...,nc-1 (indexFrom1==false) or 1,2,3,...,nc (indexFrom1==true)

    // get nnz (number of nonzero elements)
    mkIndex nnz = colPtr[nc]-colPtr[0];

    // allocate memory for the output sparse matrix
    if (!patternOnly && store2 == NULL && nnz)
        store2 = new Real[nnz];
    if (rowIdx2 == NULL && nnz)
        rowIdx2 = new mkIndex[nnz];
    if (colPtr2 == NULL)
        colPtr2 = new mkIndex[nc+1u];

    // find the length of each column after permutation
    colPtr2[0] = 0;
    mkIndex *cptr2;
    if (indexFrom1)
        cptr2 = colPtr2;
    else
        cptr2 = colPtr2 + 1;
    for (mkIndex j=0; j<nc; j++)
        cptr2[cperm[j]] = colPtr[j+1u]-colPtr[j];

    // get the right column pointers
    for (mkIndex j=0; j<nc; j++)
        colPtr2[j+1u] += colPtr2[j];

    // copy the permuted columns to the right locations
    cptr2 --;
    for (mkIndex j=0; j<nc; j++) {
        mkIndex len = colPtr[j+1u]-colPtr[j];  // column length
        mkIndex stride = colPtr[j]-colPtr[0];
        const Real *str = NULL;
        if (!patternOnly)
            str = store + stride;
        const mkIndex *ri = rowIdx + stride;
        Real *str2 = NULL;
        if (!patternOnly)
            str2 = store2 + cptr2[cperm[j]];
        mkIndex *ri2 = rowIdx2 + cptr2[cperm[j]];
        if (!patternOnly) {
            while (len--) {
                // copy elements
                *str2++ = *str++;
                // copy permuted row indices
                *ri2++ = rperm[*ri++]-indexFrom1;
            }
        }
        else {
            while (len--)
                // copy permuted row indices
                *ri2++ = rperm[*ri++]-indexFrom1;
        }
    }
}

// form a tridiagonal matrix A as a sparse matrix with three diagonals ldiag[], diag[], and udiag[]
// the lower subdiagonal, the main diagonal, and the upper subdiagonal of A are formed by
//    ldiag[0,...,nrc-2], diag[0,...,nrc-1], and udiag[0,...,nrc-2], respectively
// to be precise,
//    A(i+1,i) = ldiag[i-1] for i=1,...,nrc-1,
//    A(i,i  ) =  diag[i-1] for i=1,...,nrc, and
//    A(i,i+1) = udiag[i]   for i=1,...,nrc-1
// the output sparse matrix A is in CSC format stored in sa[], rowIdx[], colPtr[]
// if sa==NULL (rowIdx==NULL, or colPtr==NULL) is given, then the required memory will be allocated;
//    otherwise it is assumed that sufficient memory has been allocated, with the address pointed to by sa (rowIdx, colPtr, respectively)
void tridiagonal_to_spmatrix(mkIndex nrc, const Real *ldiag, const Real *diag, const Real *udiag,
                             Real *&sa, mkIndex *&rowIdx, mkIndex *&colPtr) {
    // a special case, empty input
    if (nrc == 0)
        return;

    // number of nonzero elements
    mkIndex nnz = 3u*nrc-2u;

    if (sa == NULL && nnz) {
        // allocate memory
        sa = new Real[nnz];
    }
    if (rowIdx == NULL && nnz) {
        // allocate memory
        rowIdx = new mkIndex[nnz];
    }
    if (colPtr == NULL) {
        // allocate memory
        colPtr = new mkIndex[nrc+1u];
    }

    // set the tridiagonal matrix
    Real *a0 = sa;
    mkIndex *ri = rowIdx;
    mkIndex *cp = colPtr;
    if (nrc == 1) {
        // a special case, a 1-by-1 matrix
        *a0 = *diag;
        *ri = 0;
        *cp++ = 0;
        *cp = 1;
    }
    else {
        // set A(:,1) with 2 entries
        *a0++ = *diag++;   // set value of A(1,1)
        *a0++ = *ldiag++;  // set value of A(2,1)
        *ri++ = 0;         // set row index of A(1,1)
        *ri++ = 1;         // set row index of A(2,1)
        *cp++ = 0;         // anchor of column pointers
        *cp++ = 2;         // set column pointer 1
        for (mkIndex j=1; j<nrc-1; j++) {
            // set A(:,j+1) with 3 entries
            *a0++ = *udiag++;  // set value of A(j  ,j+1)
            *a0++ = *diag++;   // set value of A(j+1,j+1)
            *a0++ = *ldiag++;  // set value of A(j+2,j+1)
            *ri++ = j-1u;      // set index of A(j  ,j+1)
            *ri++ = j;         // set index of A(j+1,j+1)
            *ri++ = j+1u;      // set index of A(j+2,j+1)
            *cp++ = 3u*j+2u;   // set column pointer j+1
        }
        // set the last column A(:,nrc) with 2 entries
        *a0++ = *udiag;  // set value of A(n-1,n)
        *a0 = *diag;     // set value of A(n-1,n)
        *ri++ = nrc-2;   // set index of A(nrc-1,nrc)
        *ri = nrc-1;     // set index of A(nrc  ,nrc)
        *cp = 3*nrc-2;   // set column pointer nrc
    }
}

#ifdef USE_MEX
// convert *mxArray to a sparse matrix
void mxArray_to_spmatrix(const mxArray *mA, mkIndex &nr, mkIndex &nc, Real *&store, mkIndex *&rowIdx, mkIndex *&colPtr) {
    if (!mxIsNumeric(mA))
        mexErrMsgTxt("mxArray_to_spmatrix(const mxArray*, mkIndex&, mkIndex&, Real&*, mkIndex*&, mkIndex*&): the input mxArray* contains non-numeric data!");
    if (!mxIsSparse(mA))
        mexErrMsgTxt("mxArray_to_spmatrix(const mxArray*, mkIndex&, mkIndex&, Real&*, mkIndex*&, mkIndex*&): the input mxArray* is not a sparse matrix!");

    // OCTAVE/MATLAB: [nr, nc] = size(A);
    nr = (mkIndex)mxGetM(mA);
    nc = (mkIndex)mxGetN(mA);

    // OCTAVE/MATLAB: [I, J, V] = find(A);
    mwIndex *ii = mxGetIr(mA);
    mwIndex *jj = mxGetJc(mA);

    // number of nonzero elements
    mkIndex nnz = jj[nc] - jj[0];

    if (rowIdx == NULL && nnz) {
        // allocate memory for row indices
        rowIdx = new mkIndex[nnz];
    }
    if (colPtr == NULL) {
        // allocate memory for column pointers
        colPtr = new mkIndex[nc+1u];
    }
    if (store == NULL && nnz) {
        // allocate memory for numerical data
        store = new Real[nnz];
    }

    // copy the row indices
    mkIndex z = nnz;
    mkIndex *ri = rowIdx;
    mwIndex *ij = ii;
    while (z--)
        *ri++ = *ij++;

    // copy the column pointers
    mkIndex c = nc+1u;
    mkIndex *cp = colPtr;
    ij = jj;
    while (c--)
        *cp++ = *ij++;

    // copy the numerical data
    if (mxIsDouble(mA)) {
        double *vv = mxGetPr(mA);
        #ifdef USE_SINGLE
            mexWarnMsgTxt("mxArray_to_spmatrix(const mxArray*, mkIndex&, mkIndex&, Real&*, mkIndex*&, mkIndex*&): the input mxArray* contains double precision floating-point data, but the output will be in single precision!");
            z = nnz;
            Real *sa = store;
            while (z--)
                *sa++ = (Real)(*vv++);
        #else
            memory_xcopy(nnz, vv, store);
        #endif
    }
    else if (mxIsSingle(mA)) {
        float *vv = (float *)mxGetData(mA);
        #ifndef USE_SINGLE
            mexWarnMsgTxt("mxArray_to_spmatrix(const mxArray*, mkIndex&, mkIndex&, Real&*, mkIndex*&, mkIndex*&): the input mxArray* contains single precision floating-point data, but the output will be in double precision!");
            z = nnz;
            Real *sa = store;
            while (z--)
                *sa++ = (Real)(*vv++);
        #else
            memory_xcopy(nnz, vv, store);
        #endif
    }
    else {
        mexErrMsgTxt("mxArray_to_spmatrix(const mxArray*, mkIndex&, mkIndex&, Real&*, mkIndex*&, mkIndex*&): the input mxArray* contains integer-type data; cannot convert it!");
    }
    if (mxIsComplex(mA))
        mexWarnMsgTxt("mxArray_to_spmatrix(const mxArray*, mkIndex&, mkIndex&, Real&*, mkIndex*&, mkIndex*&): the input mxArray* is a complex sparse matrix; the imaginary part is ignored!");
}

// convert a sparse matrix to *mxArray
mxArray *spmatrix_to_mxArray(mkSignedIndex nr, mkSignedIndex nc, const Real *store, const mkIndex *rowIdx, const mkIndex *colPtr) {
    #ifdef USE_SINGLE
        mexWarnMsgTxt("spmatrix_to_mxArray(mkSignedIndex, mkSignedIndex, const Real*, const mkIndex*, const mkIndex*): the input sparse matrix is in single precision, but the converted mxArray* will be in double precision!");
    #endif

    // create the mxArray and get its pointers
    mkIndex nnz = colPtr[nc]-colPtr[0];
    mxArray *mA = mxCreateSparse((mwSize)nr, (mwSize)nc, (mwSize)nnz, mxREAL);
    double  *sa = mxGetPr(mA);
    mwIndex *ii = mxGetIr(mA);
    mwIndex *jj = mxGetJc(mA);

    // copy the data
    while (nnz--) {
        // copy the row index
        *ii++ = (mwIndex)(*rowIdx++);
        // copy the numerical value
        *sa++ = (Real)(*store++);
    }
    nc ++;
    while (nc--) {
        // copy the column pointer
        *jj++ = (mwIndex)(*colPtr++);
    }

    return mA;
}
#endif  // of #ifdef USE_MEX



////////////////////////////////////////////////////////////////////////////
//    vector - vector operations
////////////////////////////////////////////////////////////////////////////

// concatenate two vectors v and w to form u, where
//    v, w are the input vectors of lengths nv,nw stored in sv[], sw[], respectively, and
//    u is the output vector u length nv+nw stored in su[]
// to be precise, this routine computes u(i)=v(i) for i=1,...,nv and u(nv+j)=w(j) for j=1,...,nw
// if su==NULL is given, then the required memory of size nv+nw will be allocated for the output vector u;
//    otherwise it is assumed that sufficient memory has been allocated, with the address pointed to by su
void vector_concatenate(mkIndex nv, const Real *sv, mkIndex nw, const Real *sw, Real *&su) {
    if (su == NULL && (nv || nw)) {
        // allocate memory
        su = new Real[nv+nw];
    }

    // copy elements
    memory_xcopy(nv, sv, su);
    memory_xcopy(nw, sw, su+nv);
}

// return the inner product v'*w of two vectors v and w, where
//    v, w are the input vectors of length len stored in sv[], sw[], respectively
// to be precise, the return value is v(1)*w(1)+...+v(len)*w(len)
// this routine invokes the level 1 BLAS routine xDOT if USE_BLAS is defined
Real vector_inner_product(mkIndex len, const Real *sv, const Real *sw) {
#ifdef USE_BLAS
    // invoke xDOT_
    int len2 = (int)len;
    int inc = 1;
    #ifdef USE_SINGLE
        return sdot_(&len2, sv, &inc, sw, &inc);
    #else
        return ddot_(&len2, sv, &inc, sw, &inc);
    #endif
#else
    // zero initialization
    Real ip = 0.0;

    // compute the inner product
    while (len--)
        ip += (*sv++) * (*sw++);

    // return the result
    return ip;
#endif
}

// compute A = alp*v*w' (if init0==true or sa==NULL) or A += alp*v*w' (if init0==false and sa!=NULL), where
//    alp is a scalar,
//    v, w are the input vectors of lengths nr, nc stored in sv[], sw[], respectively, and
//    A is the nr-by-nc output general matrix stored in sa[]
// if sa==NULL is given, then the required memory of size nr*nc will be allocated for the output matrix A;
//    otherwise it is assumed that sufficient memory has been allocated, with the address pointed to by sa
// this routine invokes the level 2 BLAS routine xGER if USE_BLAS is defined
void vector_outer_product(mkIndex nr, mkIndex nc, Real alp, const Real *sv, const Real *sw, Real *&sa, bool init0) {
    if (nr == 0 || nc == 0) {
        // a special case, empty input, do nothing
        return;
    }

    if (sa == NULL) {
        // allocate memory
        sa = new Real[nr*nc];
        init0 = true;
    }

#ifdef USE_BLAS
    if (init0) {
        // compute A = 0
        mkIndex sz = nr*nc;
        Real *a0 = sa;
        while (sz--)
            *a0++ = 0.0;
    }

    // invoke xGER_ for A += alp*v*w'
    int nr2 = (int)nr;
    int nc2 = (int)nc;
    int inc = 1;
    #ifdef USE_SINGLE
        sger_(&nr2, &nc2, &alp, sv, &inc, sw, &inc, sa, &nr2);
    #else
        dger_(&nr2, &nc2, &alp, sv, &inc, sw, &inc, sa, &nr2);
    #endif
#else
    if (alp != 0.0) {
        // compute A = alp*v*w' (if init0==true) or A += alp*v*w' (if init0==false)
        Real *a0 = sa;
        mkIndex jj = nc;
        while (jj--) {
            Real val = alp*(*sw++);
            mkIndex ii = nr;
            const Real *v0 = sv;
            if (init0) {
                // A(:,j)  = alp*v*w(j), where j = nc-jj
                while (ii--)
                    *a0++ = (*v0++)*val;
            }
            else {
                // A(:,j) += alp*v*w(j), where j = nc-jj
                while (ii--)
                    *a0++ += (*v0++)*val;
            }
        }
    }
    else if (init0) {
        // compute A = 0
        mkIndex sz = nr*nc;
        Real *a0 = sa;
        while (sz--)
            *a0++ = 0.0;
    }
#endif
}

// compute A = alp*v*v' (if init0==true or sa==NULL) or A += alp*v*v' (if init0==false and sa!=NULL), where
//    alp is a scalar,
//    v is the input vector of length nrc stored in sv[], and
//    A is the nrc-by-nrc output symmetric matrix stored in packed form in sa[]
// if sa==NULL is given, then the required memory of size nrc*(nrc+1)/2 for will be allocated the output symmetric matrix A;
//    otherwise it is assumed that sufficient memory has been allocated, with the address pointed to by sa
// this routine invokes the level 2 BLAS routine xSPR if USE_BLAS is defined
void vector_outer_product(mkIndex nrc, Real alp, const Real *sv, Real *&sa, bool init0) {
    if (nrc == 0) {
        // a special case, empty input, do nothing
        return;
    }

    if (sa == NULL) {
        // allocate memory
        sa = new Real[nrc*(nrc+1u)/2u];
        init0 = true;
    }

#ifdef USE_BLAS
    if (init0) {
        // zero initialization
        Real *s0 = sa;
        mkIndex sz = nrc*(nrc+1u)/2u;
        while (sz--)
            *s0++ = 0.0;
    }
    // invoke xSPR_
    char uplo = 'L';
    int inc = 1;
    int nrc2 = (int)nrc;
    #ifdef USE_SINGLE
        sspr_(&uplo, &nrc2, &alp, sv, &inc, sa);
    #else
        dspr_(&uplo, &nrc2, &alp, sv, &inc, sa);
    #endif
#else
    Real *a0 = sa;
    const Real *sw = sv;
    mkIndex jj = nrc;
    if (alp != 0.0) {
        // compute A = alp*v*v' (if init0==true) or A += alp*v*v' (if init0==false)
        while (jj) {
            mkIndex ii = jj;
            const Real *v0 = sw;
            Real val = alp*(*sw++);
            if (init0) {
                // compute A(j:nrc,j)  = alp*v(j:nrc)*v(j)
                while (ii--)
                    *a0++ = (*v0++)*val;
            }
            else {
                // compute A(j:nrc,j) += alp*v(j:nrc)*v(j)
                while (ii--)
                    *a0++ += (*v0++)*val;
            }
            jj--;
        }
    }
    else if (init0) {
        // compute A = 0
        mkIndex sz = nrc*(nrc+1u)/2u;
        while (sz--)
            *a0++ = 0.0;
    }
#endif
}

// compute A = alp*(v*w'+w*v') (if init0==true or sa==NULL) or A += alp*(v*w'+w*v') (if init0==false and sa!=NULL), where
//    alp is a scalar, and
//    v, w are the input vectors of length nrc stored in sv[], sw[], respectively
//    A is the nrc-by-nrc output symmetric matrix stored in packed form in sa[]
// if sa==NULL is given, then the required memory of size len*(len+1)/2 will be allocated for the output symmetric matrix A;
//    otherwise it is assumed that sufficient memory has been allocated, with the address pointed to by sa
// this routine invokes the level 2 BLAS routine xSPR2 if USE_BLAS is defined
void vector_symmetric_outer_product(mkIndex nrc, Real alp, const Real *sv, const Real *sw, Real *&sa, bool init0) {
    if (nrc == 0) {
        // a special case, empty input, do nothing
        return;
    }

    if (sa == NULL) {
        // allocate memory
        sa = new Real[nrc*(nrc+1u)/2u];
        init0 = true;
    }

#ifdef USE_BLAS
    if (init0) {
        // zero initialization
        mkIndex sz = nrc*(nrc+1u)/2u;
        Real *a0 = sa;
        while (sz--)
            *a0++ = 0.0;
    }
    // invoke xSPR2_
    char uplo = 'L';
    int inc = 1;
    int nrc2 = (int)nrc;
    #ifdef USE_SINGLE
        sspr2_(&uplo, &nrc2, &alp, sv, &inc, sw, &inc, sa);
    #else
        dspr2_(&uplo, &nrc2, &alp, sv, &inc, sw, &inc, sa);
    #endif
#else
    if (alp != 0.0) {
        // compute A = alp*(v*w'+w*v') (if init0==true) or A += alp*(v*w'+w*v') (if init0==false)
        Real *a0 = sa;
        mkIndex jj = nrc;
        const Real *sw2 = sw;
        const Real *sv2 = sv;
        while (jj) {
            // deal with A(j:nrc,j), where j=nrc-jj+1
            mkIndex ii = jj;
            const Real *v0 = sv2;
            const Real *w0 = sw2;
            Real vv = alp*(*sv2++);
            Real ww = alp*(*sw2++);
            if (init0) {
                // compute A(j:nrc,j) = alp*(v(j:nrc)*w(j) + w(j:nrc)*v(j))
                while (ii--)
                    *a0++  = (*v0++)*ww + (*w0++)*vv;
            }
            else {
                // compute A(j:nrc,j) = alp*(v(j:nrc)*w(j) + w(j:nrc)*v(j))
                while (ii--)
                    *a0++ += (*v0++)*ww + (*w0++)*vv;
            }
            jj--;
        }
    }
    else if (init0) {
        // compute A = 0
        mkIndex sz = nrc*(nrc+1u)/2u;
        Real *a0 = sa;
        while (sz--)
            *a0++ = 0.0;
    }
#endif
}



////////////////////////////////////////////////////////////////////////////
//    matrix - vector operations
////////////////////////////////////////////////////////////////////////////

// compute w = A*v (if init0==true or sw==NULL) or w += A*v (if init0==false and sw!=NULL), where
//    A is the nr-by-nc input general matrix stored in sa[],
//    v is the input vector of length nc stored in sv[], and
//    w is the output vector of length nr stored in sw[]
// if sw==NULL is given, then the required memory of size nr will be allocated for the output vector w;
//    otherwise it is assumed that sufficient memory has been allocated, with the address pointed to by sw
// this routine invokes the level 2 BLAS routine xGEMV if USE_BLAS is defined
void matrix_vector_multiplication(mkIndex nr, mkIndex nc, const Real *sa, const Real *sv, Real *&sw, bool init0) {
    if (nr == 0) {
        // a special case, empty output, do nothing
        return;
    }

    if (sw == NULL) {
        sw = new Real[nr];
        init0 = true;
    }
#ifdef USE_BLAS
    // level 2 BLAS routine dgemv() is about 20% faster than the homemade one in a test A*v with A of size 10000-by-10000
    char trans = 'N';
    Real alp = 1.0, bet;
    if (init0)
        bet = 0.0;
    else
        bet = 1.0;
    int inc = 1;
    int nr2 = (int)nr;
    int nc2 = (int)nc;
    #ifdef USE_SINGLE
        sgemv_(&trans, &nr2, &nc2, &alp, sa, &nr2, sv, &inc, &bet, sw, &inc);
    #else
        dgemv_(&trans, &nr2, &nc2, &alp, sa, &nr2, sv, &inc, &bet, sw, &inc);
    #endif
#else
    if (init0) {
        // compute w = 0
        Real *w0 = sw;
        mkIndex ii = nr;
        while (ii--)
            *w0++ = 0.0;
    }

    mkIndex jj = nc;
    while (jj--) {
        // compute w = w + A(:,j)*v(j), where j = nc-jj
        Real val = *sv++;
        if (val != 0.0) {
            Real *w0 = sw;
            mkIndex ii = nr;
            while (ii--) {
                // compute w(i) = w(i) + A(i,j)*v(j), where i = nr-ii
                (*w0++) += val*(*sa++);
            }
        }
        else {
            sa += nr;
        }
    }
#endif
}

// compute w' = v'*A (if init0==true or sw==NULL) or w' += v'*A (if init0==false and sw!=NULL), where
//    A is the nr-by-nc input general matrix stored in sa[],
//    v is the input vector of length nr stored in sv[], and
//    w is the output vector of length nc stored in sw[]
// if sw==NULL is given, then the required memory of size nc will be allocated for the output vector w;
//    otherwise it is assumed that sufficient memory has been allocated, with the address pointed to by sw
// this routine invokes the level 2 BLAS routine xGEMV if USE_BLAS is defined
void vector_matrix_multiplication(mkIndex nr, mkIndex nc, const Real *sv, const Real *sa, Real *&sw, bool init0) {
    if (nc == 0) {
        // a special case, empty output, do nothing
        return;
    }

    if (sw == NULL) {
        // allocate memory
        sw = new Real[nc];
        init0 = true;
    }

#ifdef USE_BLAS
    // in a test, the standard level 2 BLAS routine dgemv is about 20%-25% faster than the own implementation A*v with A of size 5000-by-5000
    char trans = 'T';
    Real alp = 1.0, bet;
    if (init0)
        bet = 0.0;
    else
        bet = 1.0;
    int inc = 1;
    int nr2 = (int)nr;
    int nc2 = (int)nc;
    #ifdef USE_SINGLE
        sgemv_(&trans, &nr2, &nc2, &alp, sa, &nr2, sv, &inc, &bet, sw, &inc);
    #else
        dgemv_(&trans, &nr2, &nc2, &alp, sa, &nr2, sv, &inc, &bet, sw, &inc);
    #endif
#else
    mkIndex jj = nc;
    Real *w0 = sw;
    while (jj--) {
        // compute w(j) = v'*A(:,j) (if init0==true) or w(j) += v'*A(:,j) (if init0==false), where j = nc-jj
        if (init0)
            *w0 = 0.0;
        mkIndex ii = nr;
        const Real *v0 = sv;
        while (ii--) {
            // compute w(j) += v(i)*A(i,j), where i = nr-ii
            *w0 += (*v0++) * (*sa++);
        }
        w0 ++;
    }
#endif
}



////////////////////////////////////////////////////////////////////////////
//    symmetric matrix - vector operations
////////////////////////////////////////////////////////////////////////////

// compute w = A*v (if init0==true or sw==NULL) or w' += v'*A (if init0==false and sw!=NULL), where
//    A is the nrc-by-nrc input symmetric matrix stored in packed form in sa[],
//    v is the input vector of length nrc stored in sv[], and
//    w is the output vector of length nrc stored in sw[]
// if sw==NULL is given, then the required memory of size nrc will be allocated for the output vector w;
//    otherwise it is assumed that sufficient memory has been allocated, with the address pointed to by sw
// this routine invokes the level 2 BLAS routine xSPMV if USE_BLAS is defined
void symmatrix_vector_multiplication(mkIndex nrc, const Real *sa, const Real *sv, Real *&sw, bool init0) {
    if (nrc == 0) {
        // a special case, empty output, do nothing
        return;
    }

    if (sw == NULL) {
        // allocate memory
        sw = new Real[nrc];
        init0 = true;
    }

#ifdef USE_BLAS
    char uplo = 'L';
    Real alp = 1.0, bet;
    if (init0)
        bet = 0.0;  // for w  = A*v
    else
        bet = 1.0;  // for w += A*v
    int inc = 1;
    int nrc2 = (int)nrc;
    #ifdef USE_SINGLE
        sspmv_(&uplo, &nrc2, &alp, sa, sv, &inc, &bet, sw, &inc);
    #else
        dspmv_(&uplo, &nrc2, &alp, sa, sv, &inc, &bet, sw, &inc);
    #endif
#else
    // in a test of multiplication of a 10000-by-10000 symmetric matrix, the own implementation is more than twice faster than the level 2 BLAS routine dspmv()
    Real *w0 = sw;
    mkIndex jj = nrc;
    if (init0) {
        // compute w = 0
        while (jj--)
            *w0++ = 0.0;
        w0 = sw;
        jj = nrc;
    }

    // compute w += A*v
    while (jj--) {
        // *sa is the j-th diagonal element, where j=nrc-jj
        *w0 += (*sv)*(*sa++);
        // now deal with off-diagonal elements of j-th row/column
        const Real *v1 = sv+1;
        Real *w1 = w0+1;
        mkIndex ii = jj;
        Real val1 = 0.0;
        Real val2 = *sv;
        while (ii--) {
            val1 += (*sa)*(*v1++);
            *(w1++) += (*sa++)*val2;
        }
        *w0 += val1;
        sv ++;
        w0 ++;
    }
#endif
}



////////////////////////////////////////////////////////////////////////////
//    sparse matrix - vector operations
////////////////////////////////////////////////////////////////////////////


// compute w = A*v (if init0==true or sw==NULL) or w += A*v (if init0==false and sw!=NULL), where
//    A is the nr-by-nc input sparse matrix in CSC format stored in store[], rowIdx[], colPtr[],
//    v is the input vector of length nc stored in sv[], and
//    w is the output vector of length nr stored in sw[]
// if sw==NULL is given, then the required memory of size nr will be allocated for the output vector w;
//    otherwise it is assumed that sufficient memory has been allocated, with the address pointed to by sw
// this routine can be used to compute w' = v'*A  or w' += v'*A, where A is a sparse matrix stored in CSR format
// this routine invokes the Intel MKL SPBLAS routine MKL_xCSCMV if USE_SPBLAS is defined
void spmatrix_vector_multiplication(mkIndex nr, mkIndex nc, const Real *store, const mkIndex *rowIdx, const mkIndex *colPtr, const Real *sv, Real *&sw, bool init0) {
    if (nr == 0) {
        // a special case, empty output, do nothing
        return;
    }
    if (sw == NULL) {
        // allocate memory
        sw = new Real[nr];
        init0 = true;
    }

#ifdef USE_SPBLAS
    char transa = 'N';
    char matdescra[] = { 'G', ' ', ' ', 'C' };  // 'G' for general, 'C' for zero-based indexing
    Real alp = 1.0, bet;
    if (init0)
        bet = 0.0;
    else
        bet = 1.0;
    int nr2 = (int)nr;
    int nc2 = (int)nc;
    #ifdef USE_SINGLE
        mkl_scscmv_(&transa, &nr2, &nc2, &alp, matdescra, store, (int *)rowIdx, (int *)colPtr, (int *)(colPtr+1), sv, &bet, sw);
    #else
        mkl_dcscmv_(&transa, &nr2, &nc2, &alp, matdescra, store, (int *)rowIdx, (int *)colPtr, (int *)(colPtr+1), sv, &bet, sw);
    #endif
#else
    if (init0) {
        // compute w = 0
        mkIndex ii = nr;
        Real *ww = sw;
        while (ii--)
            *ww++ = 0.0;
    }
    mkIndex jj = nc;
    while (jj--) {
        // compute w = w + A(:,j)*v(j), where j = nc-jj
        mkIndex kk = *(colPtr+1) - *colPtr;
        colPtr ++;
        Real val = *sv++;
        while (kk--) {
            // compute w(k) = w(k) + A(k,j)*v(j), where k = 1+(*rowIdx)
            sw[*rowIdx++] += val*(*store++);
        }
    }
#endif
}

// compute w' = v'*A (if init0==true or sw==NULL) or w' += v'*A (if init0==false and sw!=NULL), where
//    A is the nr-by-nc input sparse matrix in CSC format stored in store[], rowIdx[], colPtr[],
//    v is the input vector of length nr stored in sv[], and
//    w is the output vector of length nc stored in sw[]
// note that the number of rows in A, nr, is not required in the computation, so it is not passed
// if sw==NULL is given, then the required memory of size nc will be allocated for the output vector w;
//    otherwise it is assumed that sufficient memory has been allocated, with the address pointed to by sw
// this routine can be used to compute w = A*v  or w += A*v, where A is a sparse matrix stored in CSR format
// this routine invokes the Intel MKL SPBLAS routine MKL_xCSRMV if USE_SPBLAS is defined
void vector_spmatrix_multiplication(mkIndex nc, const Real *sv, const Real *store, const mkIndex *rowIdx, const mkIndex *colPtr, Real *&sw, bool init0) {
    if (nc == 0) {
        // a special case, empty output, do nothing
        return;
    }
    if (sw == NULL) {
        // allocate memory
        sw = new Real[nc];
        init0 = true;
    }

#ifdef USE_SPBLAS
    char transa = 'N';
    char matdescra[] = { 'G', ' ', ' ', 'C' };  // 'G' for general, 'C' for zero-based indexing
    Real alp = 1.0, bet;
    if (init0)
        bet = 0.0;
    else
        bet = 1.0;
    int nr2 = (int)nc;  // will not be used
    int nc2 = (int)nc;
    #ifdef USE_SINGLE
        mkl_scsrmv_(&transa, &nr2, &nc2, &alp, matdescra, store, (int *)rowIdx, (int *)colPtr, (int *)(colPtr+1), sv, &bet, sw);
    #else
        mkl_dcsrmv_(&transa, &nr2, &nc2, &alp, matdescra, store, (int *)rowIdx, (int *)colPtr, (int *)(colPtr+1), sv, &bet, sw);
    #endif
#else
    Real *sw2 = sw;
    mkIndex jj = nc;
    while (jj--) {
        // compute w(j) = v'*A(:,j) if init0==true, or w(j) += v'*A(:,j) if init0==false, where j = nc-jj
        mkIndex kk = *(colPtr+1) - *colPtr;
        colPtr ++;
        Real val = 0.0;
        while (kk--) {
            // compute val += v(k)*A(k,j), where k = 1+(*rowIdx)
            val += sv[*rowIdx++]*(*store++);
        }
        if (init0)
            *sw2++ = val;
        else
            *sw2++ += val;
    }
#endif
}



////////////////////////////////////////////////////////////////////////////
//    matrix - matrix operations
////////////////////////////////////////////////////////////////////////////

// horizontal concatenation of two matrices A and B, where
//    A, B are the nr-by-nca, nr-by-ncb input general matrices stored in sa[], sb[], respectively
// A and B must have the same number of rows
// the result is an nr-by-(nca+ncb) general matrix C stored in sc[]
// if sc==NULL is given, then the required memory of size nr*(nca+ncb) will be allocated for the output matrix C;
//    otherwise it is assumed that sufficient memory has been allocated, with the address pointed to by sc
// this static function resembles [A, B] or equivalently horzcat(A,B) in OCTAVE/MATLAB
void horizontal_concatenate(mkIndex nr, mkIndex nca, const Real *sa, mkIndex ncb, const Real *sb, Real *&sc) {
    if (nr == 0 || (nca == 0 && ncb == 0)) {
        // a special case, empty output, do nothing
        return;
    }

    // sizes of the input matrices
    mkIndex sza = nr*nca;
    mkIndex szb = nr*ncb;

    if (sc == NULL) {
        // allocate memory
        sc = new Real[sza+szb];
    }

    // compute C = [A, B], or equivalently C = horzcat(A,B)
    memory_xcopy(sza, sa, sc);
    memory_xcopy(szb, sb, sc+sza);
}

// vertical concatenation of two matrices A and B, where
//    A, B are the nra-by-nc, nrb-by-nc input general matrices stored in sa[], sb[], respectively
// A and B must have the same number of columns
// the result is an (nra+nrb)-by-nc general matrix C stored in sc[]
// if sc==NULL is given, then the required memory of size (nra+nrb)*nc will be allocated for the output matrix C;
//    otherwise it is assumed that sufficient memory has been allocated, with the address pointed to by sc
// this static function resembles [A; B] or equivalently vertcat(A,B) in OCTAVE/MATLAB
void vertical_concatenate(mkIndex nra, mkIndex nc, const Real *sa, mkIndex nrb, const Real *sb, Real *&sc) {
    if (nc == 0 || (nra == 0 && nrb == 0)) {
        // a special case, empty output, do nothing
        return;
    }

    if (sc == NULL) {
        // allocate memory
        sc = new Real[(nra+nrb)*nc];
    }

    // compute C = [A; B], or equivalently C = vertcat(A,B)
    Real *c0 = sc;
    while (nc--) {
        // compute C(1:nra,j) = A(:,j), where j = ncols-nc, with ncols the number of columns of A,B,C
        memory_xcopy(nra, sa, c0);
        sa += nra;
        c0 += nra;
        // compute C(1:nrb,j), where j = ncols-nc
        memory_xcopy(nrb, sb, c0);
        sb += nrb;
        c0 += nrb;
    }
}

// sum up two general matrices A,B, with zeros padded if A and B are of different dimensions, where
//    A, B are the nra-by-nca, nrb-by-ncb input general matrices stored in sa[], sb[], respectively
// the result is an nrc-by-ncc general matrix C stored in sc[], where nrc=max(nra,nrb) and ncc=max(nrb,ncb)
// if sc==NULL is given, then the required memory of size nrc*ncc will be allocated for output matrix C;
//    otherwise, it is assumed that sufficient memory has been allocated, with the address pointed to by sc
void matrix_xsum(mkIndex nra, mkIndex nca, const Real *sa, mkIndex nrb, mkIndex ncb, const Real *sb, Real *&sc) {
    // find dimensions
    mkIndex MM = (nra>=nrb) ? nra : nrb;
    mkIndex NN = (nca>=ncb) ? nca : ncb;

    if (MM == 0 || NN == 0) {
        // a special case, empty output, do nothing
        return;
    }

    if (sc == NULL) {
        // allocate memory
        sc = new Real[MM*NN];
    }

    // min of numbers of rows and min of numbers of columns
    mkIndex mm = (nra<=nrb) ? nra : nrb;
    mkIndex nn = (nca<=ncb) ? nca : ncb;

    // compute C(:,1:nn) = A(:,1:nn) + B(:,1:nn)
    Real *c0 = sc;
    mkIndex jj = nn;
    while (jj--) {
        // compute C(:,j) = A(:,j) + B(:,j), where j = nn-jj
        mkIndex ii = mm;
        while (ii--)
            (*c0++) = (*sa++) + (*sb++);
        if (nra > mm) {
            // A has more rows than B; compute C(mm+1:nra,j) = A(mm+1:nra,j)
            ii = nra-mm;
            while (ii--)
                *c0++ = *sa++;
        }
        else if (nrb > mm) {
            // B has more rows than A; compute C(mm+1:nrb,j) = B(mm+1:nrb,j)
            ii = nrb-mm;
            while (ii--)
                *c0++ = *sb++;
        }
    }

    // dealing with the last columns, if any
    if (nca > nn) {
        // A has more columns than B; compute C(:,nca+1:nn) = A(:,nca+1:nn)
        mkIndex jj = nca - nn;
        while (jj--) {
            // compute C(:,j) = A(:,j), where j = nca-jj
            mkIndex ii = nra;
            while (ii--) {
                // compute C(1:nra,j) = A(1:nra,j)
                *c0++ = *sa++;
            }
            if (nra < MM) {
                // C has more rows than A, padding zeros
                ii = MM - nra;
                while (ii--)
                    *c0++ = 0.0;
            }
        }
    }
    else if (ncb > nn) {
        // B has more columns than A; compute C(:,ncb+1:nn) = B(:,ncb+1:nn)
        mkIndex jj = ncb - nn;
        while (jj--) {
            // compute C(:,j) = B(:,j), where j = ncb-jj
            mkIndex ii = nrb;
            while (ii--) {
                // compute C(1:nra,j) = A(1:nra,j)
                *c0++ = *sb++;
            }
            if (nrb < MM) {
                // C has more rows than B, padding zeros
                ii = MM - nrb;
                while (ii--)
                    *c0++ = 0.0;
            }
        }
    }
}

// compute C = A*B (if init0==true or sc==NULL) or C += A*B (if init0==false and sc!=NULL), where
//    A, B are the nr-by-nrc, nrc-by-nc input general matrices stored in sa[], sb[], respectively, and
//    C is the nr-by-nc output general matrix stored in sc[]
// if sc==NULL is given, then the required memory nr*nc will be allocated for the output matrix C;
//    otherwise it is assumed that sufficient memory has been allocated, with the address pointed to by sc
// this routine invokes the level 3 BLAS routine xGEMM if USE_BLAS is defined
void matrix_matrix_multiplication(mkIndex nr, mkIndex nrc, mkIndex nc, const Real *sa, const Real *sb, Real *&sc, bool init0) {
    if (nr == 0 || nc == 0) {
        // a special case, empty A or B, do nothing
        return;
    }

    if (sc == NULL) {
        // allocate memory
        sc = new Real[nr*nc];
        init0 = true;
    }

#ifdef USE_BLAS
    // in a test of multiplication of two matrices of size 5000-by-5000, level 3 BLAS routine dgemm is about 6-7 times faster
    char trans = 'N';
    Real alp = 1.0, bet;
    if (init0)
        bet = 0.0;
    else
        bet = 1.0;
    int m = (int)nr;
    int k = (int)nrc;
    int n = (int)nc;
    #ifdef USE_SINGLE
        sgemm_(&trans, &trans, &m, &n, &k, &alp, sa, &m, sb, &k, &bet, sc, &m);
    #else
        dgemm_(&trans, &trans, &m, &n, &k, &alp, sa, &m, sb, &k, &bet, sc, &m);
    #endif
#else
    Real *sc2 = sc;
    if (init0) {
        // compute C = 0
        mkIndex sz = nr*nc;
        while (sz--)
            *sc2++ = 0.0;
        sc2 = sc;
    }
    mkIndex jj = nc;
    while (jj--) {
        // compute C(:,j) = A*B(:,j), where A is *this and j = nc-jj
        const Real *a0 = sa;
        mkIndex kk = nrc;
        while (kk--) {
            // compute C(:,j) = C(:,j) + A(:,k)*B(k,j), where k = nrc-kk
            Real val = *sb++;
            if (val != 0.0) {
                Real *c0 = sc2;
                mkIndex ii = nr;
                while (ii--) {
                    // compute C(i,j) = C(i,j) + A(i,k)*B(k,j), where i = nr-ii
                    (*c0++) += val*(*a0++);
                }
            }
            else {
                a0 += nr;
            }
        }
        sc2 += nr;
    }
#endif
}



////////////////////////////////////////////////////////////////////////////
//    matrix - symmetric matrix operations
////////////////////////////////////////////////////////////////////////////

// compute C = A + bet*B (if sa!=sc) or A += bet*B (if sa==sc), where
//    A is the nrc-by-nrc input general matrix stored in sa[],
//    B is the nrc-by-nrc input symmetric matrix stored in packed form in sb[], and
//    C is the nrc-by-nrc output general matrices stored in sc[]
// if sc==NULL is given, then the required memory of size nrc*nrc will be allocated for the output matrix C;
//    otherwise it is assumed that sufficient memory has been allocated, with the address pointed to by sc
void matrix_symmatrix_addition(mkIndex nrc, const Real *sa, Real bet, const Real *sb, Real *&sc) {
    if (nrc == 0) {
        // a special case, empty A and B, do nothing
        return;
    }

    if (sc == NULL) {
        // allocate memory
        sc = new Real[nrc*nrc];
    }

    if (sc != sa) {
        // compute C = A
        memory_xcopy(nrc*nrc, sa, sc);
    }

    if (bet != 0.0) {
        // compute C = C + bet*B
        Real *c0 = sc;
        const Real *b0 = sb;
        for (mkIndex jj=0; jj<nrc; jj++) {
            // for address of C(j,j), where j = jj+1
            c0 += jj;
            Real *c1 = c0;
            // compute C(j,j) += bet*B(j,j)
            *(c0++) += bet*(*b0++);
            // c0 points to C(j+1,j) now
            mkIndex ii = nrc-jj-1;
            // add the off-diagonal elements in column j+1 of B and the elements in column j+1 in the the strictly lower triangular part of A
            while (ii--) {
                Real val = bet*(*b0++);
                // compute C(j+i,j) += bet*B(j+i,j), where i = (nrows-jj)-ii
                *(c0++) += val;
                // compute C(j,j+i) += bet*B(j,j+i); note that B(j,j+i) = =B(j+i,j)
                *(c1+=nrc) += val;
            }
        }
    }
}

// compute C = A*B (if init0==true or sc==NULL) or C += A*B (if init0==false and sc!=NULL), where
//    A is the nr-by-nrc input general matrix stored in sa[],
//    B is the nrc-by-nrc input symmetric matrix stored in packed form in sb[], and
//    C is the nr-by-nrc output general matrices stored in sc[]
// if sc==NULL is given, then the required memory will be allocated for the output matrix C;
//    otherwise it is assumed that sufficient memory has been allocated, with the address pointed to by sc
void matrix_symmatrix_multiplication(mkIndex nr, mkIndex nrc, const Real *sa, const Real *sb, Real *&sc, bool init0) {
    if (nr == 0 || nrc == 0) {
        // a special case, empty A or B, do nothing
        return;
    }

    if (sc == NULL) {
        // allocate memory
        sc = new Real[nr*nrc];
        init0 = true;
    }

    if (init0) {
        // compute C = 0
        Real *c0 = sc;
        mkIndex sz = nr*nrc;
        while (sz--)
            *c0++ = 0.0;
    }

    // compute C += A*B
    // initial pointer to A(1,1)
    const Real *a0 = sa;
    // initial pointer to B(1,1)
    const Real *b0 = sb;
    for (mkIndex k=1; k<=nrc; k++) {
        // compute C += A(:,k)*B(k,:)
        // now a0 points to A(1,k) and b0 points to B(k,1)

        // for pointer to B(k,1)
        const Real *b1 = b0;
        // for pointer to C(1,1)
        Real *c0 = sc;
        for (mkIndex j=1; j<=nrc; j++) {
            // pointer for A(1,k)
            const Real *a2 = a0;
            // value of B(k,j)
            Real val = *b1;

            // compute C(:,j) += A(:,k)*B(k,j)
            mkIndex ii = nr;
            while (ii--) {
                // compute C(i,j) += A(i,k)*B(k,j), where i = nr-ii
                (*c0++) += (*a2++)*val;
            }

            if (j < k) {
                // pointer to B(k,j) => pointer to B(k,j+1)
                b1 += (nrc-j);
            }
            else {
                // pointer to B(j,k) => pointer to B(j+1,k)
                b1 ++;
            }
        }
        // pointer to A(1,k) => pointer to A(1,k+1)
        a0 += nr;
        // pointer to B(k,1) => pointer to B(k+1,1)
        b0 ++;
    }
}

// compute C = A*B (if init0==true or sc==NULL) or C += A*B (if init0==false and sc!=NULL), where
//    A is the nrc-by-nrc input symmetric matrix stored in packed form in sa[],
//    B is the nrc-by-nc input general matrix stored in sb[], and
//    C is the nrc-by-nc output general matrix stored in sb[]
// if sc==NULL is given, then the required memory of size nrc*nc will be allocated for the output matrix C;
//    otherwise it is assumed that sufficient memory has been allocated, with the address pointed to by sc
void symmatrix_matrix_multiplication(mkIndex nrc, mkIndex nc, const Real *sa, const Real *sb, Real *&sc, bool init0) {
    if (nrc == 0 || nc == 0) {
        // a special case, empty A or B, do nothing
        return;
    }

    if (sc == NULL) {
        // allocate memory
        sc = new Real[nrc*nc];
        init0 = true;
    }

    // pointer to C(1,1)
    Real *c0 = sc;

    if (init0) {
        // compute C = 0
        mkIndex sz = nrc*nc;
        while (sz--)
            *c0++ = 0.0;
        c0 = sc;
    }

    while (nc--) {
        // compute C(:,k) += A*B(:,k), where k = ncols-nc, with ncols the number of columns in C

        // pointer to A(1,1)
        const Real *a0 = sa;

        mkIndex jj = nrc;
        while (jj--) {
            // compute C(j,k) += A(j,:)*B(:,k) and C(j,j+1:nrc) += A(j,j+1:nrc)*B(j,k), where j = nrc-jj
            // now a0 points to A(j,j), sb points to B(j,k), and c0 points to C(j,k)

            // value of B(j,k)
            Real bjk = *sb;

            // cjk will be used for C(j,k) += cjk
            // compute cjk = A(j,j)*B(j,k)
            Real cjk = bjk*(*a0++);

            // now deal with off-diagonal elements of j-th row/column of A
            const Real *b1 = sb + 1;
            Real *c1 = c0 + 1;
            mkIndex ii = jj;
            while (ii--) {
                // compute cjk += A(j,i)*B(i,k), where i = j+(jj-ii)
                cjk += (*a0)*(*b1++);
                // compute C(i,k) += A(i,j)*B(j,k)
                *(c1++) += (*a0++)*bjk;
            }

            // pointer to B(j,k) => pointer to B(j+1,k), except for the last iteration in the jj-loop
            sb ++;

            // compute C(j,k) += cjk
            // pointer to C(j,k) => pointer to C(j+1,k), except for the last iteration in the jj-loop
            *(c0++) += cjk;
        }
    }
}



////////////////////////////////////////////////////////////////////////////
//    matrix - sparse matrix operations
////////////////////////////////////////////////////////////////////////////

// compute C = A + bet*B (if sa!=sc) or A += bet*B (if sa==sc), where
//    A is the nr-by-nc input general matrix stored in sa[],
//    B is the nr-by-nc input sparse matrix in CSC format stored in sb[], rowIdx[], colPtr[], and
//    C is the nr-by-nc output general matrices stored in sc[]
// if sc==NULL is given, then the required memory of size nr*nc will be allocated for the output matrix C;
// otherwise it is assumed that sufficient memory has been allocated, with the address pointed to by sc
void matrix_spmatrix_addition(mkIndex nr, mkIndex nc, const Real *sa, Real bet, const Real *sb, const mkIndex *rowIdx, const mkIndex *colPtr, Real *&sc) {
    if (nr == 0 || nc == 0) {
        // a special case, empty A and B, do nothing
        return;
    }

    if (sc == NULL) {
        // allocate memory
        sc = new Real[nr*nc];
    }

    if (sc != sa) {
        // compute C = A
        memory_xcopy(nr*nc, sa, sc);
    }

    // compute C = C + bet*B
    if (bet != 0.0) {
        Real *c0 = sc;
        mkIndex jj = nc;
        while (jj--) {
            // now c0 points to C(1,j), where j = nc-jj
            mkIndex kk = *(colPtr+1) - (*colPtr);
            colPtr ++;
            // compute C(:,j) = C(:,j) + bet*B(:,j)
            while (kk--) {
                // compute C(i,j) = C(i,j) + bet*B(i,j), where i=1+(*rowIdx)
                c0[*rowIdx++] += bet*(*sb++);
            }
            c0 += nc;
        }
    }
}

// compute C = A*B (if init0==true or sc==NULL) or C += A*B (if init0==false and sc!=NULL), where
//    A is the nr-by-nrc general matrix stored in sa[],
//    B is the nrc-by-nc input sparse matrix in CSC format stored in sb[], rowIdx[], colPtr[], and
//    C is the nr-by-nc output general matrix stored in sc[]
// if sc==NULL is given, then the required memory of size nr*nc will be allocated for the output matrix C;
//    otherwise it is assumed that sufficient memory has been allocated, with the address pointed to by sc
void matrix_spmatrix_multiplication(mkIndex nr, mkIndex nrc, mkIndex nc, const Real *sa,
                                    const Real *sb, const mkIndex *rowIdx, const mkIndex *colPtr,
                                    Real *&sc, bool init0) {
    if (nr == 0 || nc == 0) {
        // a special case, empty A or B, do nothing
        return;
    }

    if (sc == NULL) {
        // allocate memory
        sc = new Real[nr*nc];
        init0 = true;
    }

    // pointer to C(1,1)
    Real *c0 = sc;

    if (init0) {
        // compute C = 0
        mkIndex sz = nr*nc;
        while (sz--)
            *c0++ = 0.0;
        c0 = sc;
    }

    mkIndex jj = nc;
    while (jj--) {
        // compute C(:,j) = A*B(:,j), where j = nc-jj
        mkIndex kk = *(colPtr+1) - *colPtr;
        colPtr ++;
        while (kk--) {
            // store val = B(*rowIdx,j);
            Real val = *sb++;
            // find a0 pointing to A(1,r), where r=1+(*rowIdx)
            const Real *a0 = sa + (*rowIdx++)*nr;
            Real ii = nr;
            Real *ci = c0;
            // compute C(:,j) = C(:,j) + A(:,1+r)*B(r,j)
            while (ii--)
                *(ci++) += val*(*a0++);
        }
        c0 += nr;
    }
}

// compute C = A*B (if init0==true or sc==NULL) or C += A*B (if init0==false and sc!=NULL), where
//    A is the nr-by-nrc input sparse matrix in CSC format stored in sa[], rowIdx[], colPtr[],
//    B is the nrc-by-nc input general matrix stored in sb[], and
//    C is the nr-by-nc output general matrix stored in sc[]
// if sc==NULL is given, then the required memory of size nr*nc will be allocated for the output matrix C;
//    otherwise, it is assumed that sufficient memory has been allocated, with the address pointed to by sc
void spmatrix_matrix_multiplication(mkIndex nr, mkIndex nrc, mkIndex nc, const Real *sa, const mkIndex *rowIdx, const mkIndex *colPtr,
                                    const Real *sb, Real *&sc, bool init0) {
    if (nr == 0 || nc == 0) {
        // a special case, empty A or B, do nothing
        return;
    }

    if (sc == NULL) {
        // allocate memory for C
        sc = new Real[nr*nc];
        init0 = true;
    }

    // pointer to C(1,1)
    Real *c0 = sc;

    if (init0) {
        // compute C = 0
        mkIndex sz = nr*nc;
        while (sz--)
            *c0++ = 0.0;
        c0 = sc;
    }

    mkIndex jj = nc;
    while (jj--) {
        // compute C(:,j) += A * B(:,j), where j = nc-jj
        const Real *a0 = sa;
        const mkIndex *ri = rowIdx;
        const mkIndex *cp = colPtr;
        mkIndex ii = nrc;  // number of columns in A, as well as number of rows in B
        while (ii--) {
            // compute C(:,j) += A(:,i) * B(i,j), where i = nrc-ii
            mkIndex kk = *(cp+1) - *cp;
            cp ++;
            Real val = *sb++;
            while (kk--) {
                // compute C(k,j) += A(k,i) * B(i,j), where k = 1+(*ri)
                c0[*ri++] += val*(*a0++);
            }
        }
        c0 += nr;
    }
}



////////////////////////////////////////////////////////////////////////////
//    symmetric matrix - symmetric matrix operations
////////////////////////////////////////////////////////////////////////////

// compute C = A*B (if init0==true or sc==NULL) or C += A*B (if init0==false and sc!=NULL), where
//    A, B are the nrc-by-nrc input symmetric matrices stored in packed form in sa[] and sb[], respectively, and
//    C is the nrc-by-nrc output general matrix stored in sc[]
// if sc==NULL is given, then the required memory of size nrc*nrc will be allocated for the output matrix C;
//    otherwise it is assumed that sufficient memory has been allocated, with the address pointed to by sc
void symmatrix_symmatrix_multiplication(mkIndex nrc, const Real *sa, const Real *sb, Real *&sc, bool init0) {
    if (nrc == 0) {
        // a special case, empty A and B, do nothing
        return;
    }

    if (sc == NULL) {
        // allocate memory
        sc = new Real[nrc*nrc];
        init0 = true;
    }

    if (init0) {
        // compute C = 0
        Real *c0 = sc;
        mkIndex sz = nrc*nrc;
        while (sz--)
            *c0++ = 0.0;
    }

    // compute C += A*B
    // initial pointer to A(1,1)
    const Real *a0 = sa;
    // initial pointer to B(1,1)
    const Real *b0 = sb;
    for (mkIndex k=1; k<=nrc; k++) {
        // compute C += A(:,k)*B(k,:)
        // now a0 points to A(k,1) and b0 points to B(k,1)

        // for pointer to A(k,1)
        const Real *a1 = a0;
        // for pointer to B(k,1)
        const Real *b1 = b0;
        // for pointer to C(1,1)
        Real *c0 = sc;
        for (mkIndex j=1; j<=nrc; j++) {
            // pointer for A(k,1)
            const Real *a2 = a1;
            // value of B(k,j)
            Real val = *b1;

            // compute C(:,j) += A(:,k)*B(k,j) in 2 installments
            // the first installment, C(1:k,j) += A(1:k,k)*B(k,j)
            for (mkIndex i=1; i<k; i++) {
                // C(i,j) += A(k,i)*B(k,j), where A(k,i) == A(i,k)
                (*c0++) += (*a2)*val;
                a2 += (nrc-i);
            }
            // the second installment, C(k:nrc,j) += A(k:nrc,k)*B(k,j)
            for (mkIndex i=k; i<=nrc; i++) {
                // C(i,j) += A(i,k)*B(k,j)
                (*c0++) += (*a2++)*val;
            }

            if (j < k) {
                // for pointer to B(k,j+1)
                b1 += (nrc-j);
            }
            else {
                // for pointer to B(j+1,k)
                b1 ++;
            }
        }
        // for pointer to A(k+1,1)
        a0 ++;
        // for pointer to B(k+1,1)
        b0 ++;
    }
}



////////////////////////////////////////////////////////////////////////////
//    sparse matrix - sparse matrix operations
////////////////////////////////////////////////////////////////////////////

// return true if A == B mathematically, and otherwise return false, where
//    A is the nr-by-nc input sparse matrix in CSC format stored in sa[], rowIdx0[], colPtr0[]
//    B is the nr-by-nc input sparse matrix in CSC format stored in sb[], rowIdx1[], colPtr1[]
// the number of rows, nr, is not required in the computation, so that it is not passed
// if patternOnly==true, then only the sparsity pattern is considered, and sa[], sb[] will not be used
// this routine requires that the elements in each column be sorted with respect to the row indices for both A and B;
bool is_spmatrix_equal_spmatrix(mkIndex nc, const Real *sa, const mkIndex *rowIdx0, const mkIndex *colPtr0,
                                            const Real *sb, const mkIndex *rowIdx1, const mkIndex *colPtr1, bool patternOnly) {
    if (patternOnly) {
        mkIndex nnz = colPtr0[nc] - colPtr0[0];  // number of stored entries in A
        if (*colPtr0 == *colPtr1) {
            if (are_elements_equal_elements(nc+1u, colPtr0, colPtr1) == false)
                return false;
        }
        else {
            for (mkIndex i=1; i<=nc; i++) {
                if (colPtr0[i]-colPtr0[0] != colPtr1[i]-colPtr1[0])
                    return false;
            }
        }
        return are_elements_equal_elements(nnz, rowIdx0, rowIdx1);
    }

    // patternOnly == false
    while (nc--) {
        // compare A(:,j) with B(:,j), where j = ncols-nc, with ncols the number of columns
        mkIndex k0 = *(colPtr0+1) - *colPtr0;  // number of stored entries in the current column of A
        mkIndex k1 = *(colPtr1+1) - *colPtr1;  // number of stored entries in the current column of B
        while (k0 && k1) {
            if (*rowIdx0 == *rowIdx1) {
                 if (*sa++ != *sb++)
                     return false;
                rowIdx0 ++;
                rowIdx1 ++;
                k0 --;
                k1 --;
            }
            else if (*rowIdx0 < *rowIdx1) {
                if (*sa++ != 0.0)
                    return false;
                rowIdx0 ++;
                k0 --;
            }
            else {  // *rowIdx0 > *rowIdx1
                if (*sb++ != 0.0)
                    return false;
                rowIdx1 ++;
                k1 --;
            }
        }
        // the rest in column j (=ncols-nc, with ncols the number of columns), if any, must be zero for returning true
        while (k0--) {
            if (*sa++ != 0.0)
                return false;
            rowIdx0 ++;
        }
        while (k1--) {
            if (*sb++ != 0.0)
                return false;
            rowIdx1 ++;
        }
        // column j (=ncols-nc, with ncols the number of columns) is done; continue with the next column
        colPtr0 ++;
        colPtr1 ++;
    }
    return true;
}

// compute C = alp*A + bet*B, where
//    A is the nr-by-nc input sparse matrix in CSC format stored in sa[], rowIdx0[], colPtr0[], and
//    B is the nr-by-nc input sparse matrix in CSC format stored in sb[], rowIdx1[], colPtr1[], and
//    C is the nr-by-nc output sparse matrix in CSC format stored in sc[], rowIdx2[], colPtr2[]
// the number of rows, nr, is not required in the computation, so that it is not passed
// if sc==0 (rowIdx2==NULL, or colPtr2==NULL) is given, then sufficient memory will be allocated;
//    otherwise it is assumed that sufficient memory has been allocated, with the address pointed to by sc (rowIdx2, or colPtr2, respectively)
// if patternOnly==true, then only the sparsity pattern is considered, and alp, bet, sa[], sb[], and sc[] will not be used
// this routine requires that the elements in each column be sorted with respect to the row indices for both A and B;
//    the output matrix C will also have elements in each column being sorted
void spmatrix_spmatrix_addition(mkIndex nc,
                                Real alp, const Real *sa, const mkIndex *rowIdx0, const mkIndex *colPtr0,
                                Real bet, const Real *sb, const mkIndex *rowIdx1, const mkIndex *colPtr1,
                                Real *&sc, mkIndex *&rowIdx2, mkIndex *&colPtr2, bool patternOnly) {
    // need nc+1 column pointers in all cases
    if (colPtr2 == NULL)
        colPtr2 = new mkIndex[nc+1u];

    // deal with special cases alp==0.0 or bet==0.0 first (when patternOnly==false)
    if (!patternOnly) {
        if (alp == 0.0 && bet == 0.0) {
            // output is a zero matrix with nc columns
            for (mkIndex j=0; j<=nc; j++)
                colPtr2[j] = 0;
            return;
        }
        if (alp == 0.0) {
            // output is bet*B
            memcpy(colPtr2, colPtr1, (nc+1)*sizeof(mkIndex));
            mkIndex nnz = colPtr1[nc] - colPtr1[0];
            if (rowIdx2 == NULL && nnz)
                rowIdx2 = new mkIndex[nnz];
            memcpy(rowIdx2, rowIdx1, nnz*sizeof(mkIndex));
            if (sc == NULL && nnz)
                sc = new Real[nnz];
            if (bet == 1.0) {
                memory_xcopy(nnz, sb, sc);
            }
            else {
                Real *c0 = sc;
                while (nnz--)
                    *c0++ = bet*(*sb++);
            }
            return;
        }
        if (bet == 0.0) {
            // output is alp*A
            memcpy(colPtr2, colPtr0, (nc+1)*sizeof(mkIndex));
            mkIndex nnz = colPtr0[nc] - colPtr0[0];
            if (rowIdx2 == NULL && nnz)
                rowIdx2 = new mkIndex[nnz];
            memcpy(rowIdx2, rowIdx0, nnz*sizeof(mkIndex));
            if (sc == NULL && nnz)
                sc = new Real[nnz];
            if (alp == 1.0) {
                memory_xcopy(nnz, sa, sc);
            }
            else {
                Real *c0 = sc;
                while (nnz--)
                    *c0++ = alp*(*sa++);
            }
            return;
        }
    }

    // allocate sufficient memory
    mkIndex nnz0 = colPtr0[nc] - colPtr0[0];
    mkIndex nnz1 = colPtr1[nc] - colPtr1[0];
    mkIndex maxnnz2 = nnz0 + nnz1;
    if (!patternOnly && sc == NULL && maxnnz2)
        sc = new Real[maxnnz2];
    Real *c0 = sc;
    if (rowIdx2 == NULL && maxnnz2)
        rowIdx2 = new mkIndex[maxnnz2];

    // compute C = alp*A + bet*B
    mkIndex *ridx2 = rowIdx2;
    mkIndex *cptr2 = colPtr2;
    *cptr2 = 0;  // the anchor
    while (nc--) {
        // compute C(:,j) = alp*A(:,j) + bet*B(:,j), where j = ncols-nc, with ncols the number of columns
        mkIndex k0 = *(colPtr0+1) - *colPtr0;  // number of stored entries in the current column of A
        colPtr0 ++;
        mkIndex k1 = *(colPtr1+1) - *colPtr1;  // number of stored entries in the current column of B
        colPtr1 ++;
        mkIndex k2 = 0;
        while (k0 && k1) {
            if (*rowIdx0 == *rowIdx1) {
                if (patternOnly) {
                    *ridx2++ = *rowIdx0;
                    k2 ++;
                }
                else {
                    *c0 = alp*(*sa++) + bet*(*sb++);
                    if (*c0 != 0.0) {
                        *ridx2++ = *rowIdx0;
                        c0 ++;
                        k2 ++;
                    }
                    // else the current entry is zero; do not store it in C
                }
                rowIdx0 ++;
                rowIdx1 ++;
                k0 --;
                k1 --;
            }
            else if (*rowIdx0 < *rowIdx1) {
                if (!patternOnly)
                    *c0++ = alp*(*sa++);
                *ridx2++ = *rowIdx0++;
                k0 --;
                k2 ++;
            }
            else {  // *rowIdx0 > *rowIdx1
                if (!patternOnly)
                    *c0++ = bet*(*sb++);
                *ridx2++ = *rowIdx1++;
                k1 --;
                k2 ++;
            }
        }
        if (k0) {
            k2 += k0;
            while (k0--) {
                if (!patternOnly)
                    *c0++ = alp*(*sa++);
                *ridx2++ = *rowIdx0++;
            }
        }
        else if (k1) {
            k2 += k1;
            while (k1--) {
                if (!patternOnly)
                    *c0++ = bet*(*sb++);
                *ridx2++ = *rowIdx1++;
            }
        }
        *(cptr2+1) = *(cptr2) + k2;
        cptr2 ++;
        // column j (=ncols-nc, with ncols the number of columns) of C is done; continue with the next column
    }
}

#ifdef USE_NAMESPACE
}  // end of namespace NEWMAT
#endif
