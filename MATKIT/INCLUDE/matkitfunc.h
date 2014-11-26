#ifndef MATKIT_FUNC_H
#define MATKIT_FUNC_H
/*! \file matkitfunc.h
 *  A number of <b>MATKIT</b> functions without classes are declared here.
 *
 *  The output vectors, general matrices, symmetric matrices, and sparse matrices
 *  are generally passed by references to pointers.
 *  If a reference to a pointer pointing to <tt>NULL</tt>,
 *  then the required memory will be allocated by C++ operator <b>new</b>.
 *  Therefore, the allocated memory should be freed by the operator <b>delete</b>.
 *
 *  One can allocate memory by her or his own, e.g. the C routine <b>malloc</b>(),
 *  and pass the memory addresses to a <b>MATKIT</b> function (without classes).
 *  In this way these routines will not allocate memory and the memory management
 *  is handled outside the <b>MATKIT</b> routines (without classes).
 *
 *  Exceptions are handled, if any, by returning an <b>int</b> flag which signifies the status.
 *  For example, matrix_market_spmatrix_read() encounters a file read error,
 *  it returns a negative number.
 *
 *  Again <b>MATKIT</b> uses C++ operators <b>new</b> and <b>delete</b> for memory management.
 *  The <b>new</b> and <b>delete</b> operators can be overridden if it is preferred
 *  to use another memory management scheme.
 */

#include <stddef.h>  // for NULL pointer, pointer subtraction, etc.
#include "matkitdef.h"


#ifdef USE_NAMESPACE
namespace MATKIT {
#endif


////////////////////////////////////////////////////////////////////////////
//    miscellaneous routines
////////////////////////////////////////////////////////////////////////////

/** @defgroup group_misc Miscellaneous Routines
 *
 *  This module provides some miscellaneous routines (declared in matkitfunc.h).
 *
 *  @{
 */

//! Return <b>true</b> if <em>sa</em>[<em>i</em>]==<em>sb</em>[<em>i</em>] for <em>i=0,...,len-1</em>, and otherwise return <b>false</b>.
/*! Conceptually, this routine tells whether <em>v==w</em>, where
 *     <em>v, w</em> are the input vectors of length <em>len</em> stored in <em>sv</em>[], <em>sw</em>[], respectively.
 *  \remark  It can be used to tell whether <em>A==B</em>, where <em>A</em> and <em>B</em> are both general or both symmetric matrices.
 */
// return true if sa[i]==sb[i] for i=0,...,len-1, and otherwise return false
// conceptually, this routine tells whether v==w, where
//    v, w are the input vectors of length len stored in sv[], sw[], respectively
// it can be used to tell whether A==B, where A and B are both general or both symmetric matrices
template<class Index, class Element>
bool are_elements_equal_elements(Index len, const Element *sa, const Element *sb) {
     if (len <= 0) {
          // a special case, empty input
          return true;
     }

     while (len--) {
          if (*sa++ != *sb++)
               return false;
     }
     return true;
}

//! Compute <em>des0</em>[<em>j*inc2</em>] = <em>src0</em>[<em>j*inc1</em>] for <em>j=0,...,len-1</em>, where <em>src0=src</em> if <em>inc1>=0</em>, or <em>src0=src</em>-(<em>len-1</em>)*<em>inc1</em> if <em>inc1<0</em>, and <em>des0=des</em> if <em>inc2>=0</em>, or <em>des0=des</em>-(<em>len-1</em>)*<em>inc2</em> if <em>inc2<0</em>.
/** Conceptually, this routine computes <em>w=v</em>.
 *  <ul>
 *  <li>  <em>v</em> is the input  vector stored in <em>src0</em>[<em>0,inc1,...,(len-1)*inc1</em>].
 *  <li>  <em>w</em> is the output vector stored in <em>des0</em>[<em>0,inc2,...,(len-1)*inc2</em>].
 *  </ul>
 *
 *  \param  len     (as input) is the length of <em>v, w</em>.
 *  \param  src[]   (as input) stores the elements of <em>v</em>.
 *  \param  des[]   (as output) stores the elements of <em>w</em>.
 *  \param  inc1    (as input) is the index increment of <em>src</em>[] elements.
 *  \param  inc2    (as input) is the index increment of <em>des</em>[] elements.
 *
 *  It can be used to compute <em>A=B</em>, where <em>A</em> and <em>B</em> are both general or both symmetric matrices.
 *
 *  It can be used to copy nonzero elements of a sparse matrix.
 *
 *  It can also be used to compute a submatrix or a subvector (i.e. copy elements).
 *
 *  \remark  This routine invokes the level 1 <b>BLAS</b> routine <b>xCOPY</b> if <tt>USE_BLAS</tt> is defined.
 */
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
void memory_xcopy(mkIndex len, const Real *src, Real *des, mkSignedIndex inc1=1, mkSignedIndex inc2=1);

//! Compute <em>sv</em>[<em>i*inc</em>]<em>=alp*sv</em>[<em>i*inc</em>] for <em>i=0,...,len-1</em>.
/*! Conceptually, it computes <em>v*=alp</em>, where <em>v</em> is stored in <em>sv</em>[<em>0,inc,...,(len-1)*inc</em>].
 *
 *  \param  len     (as input) is the length of <em>v</em>.
 *  \param  alp     (as input) is a scalar.
 *  \param  sv[]    (as output) stores the elements of <em>v</em>.
 *  \param  inc     (as input) is the index increment of <em>sv</em>[] elements.
 *
 *  It can also be used to compute <em>A*=alp</em>, where <em>A</em> is a general, symmetric, or sparse matrix.
 *
 *  \remark  This routine invokes the level 1 <b>BLAS</b> routine <b>xSCAL</b> if <tt>USE_BLAS</tt> is defined.
 */
// compute sv[i*inc] *= alp for i=0,...,len-1
// conceptually, this routine computes v *= alp, where v is stored in sv[0,inc,...,(len-1)*inc]
// it can also be used to compute A *= alp, where A is a general, symmetric, or sparse matrix
// this routine invokes the level 1 BLAS routine xSCAL if USE_BLAS is defined
void in_place_scale_elements(mkIndex len, Real alp, Real *sv, mkIndex inc=1);

//! Compute <em>sw</em>[<em>i</em>]=<em>alp*sv</em>[<em>i</em>] (if <em>init0==true</em> or <em>sw==NULL</em>) or <em>sw</em>[<em>i</em>]+=<em>alp*sv</em>[<em>i</em>] (if <em>init0</em>==<b>false</b> and <em>sw</em>\f$\neq\f$<tt>NULL</tt>) for <em>i=0,...,len-1</em>.
/*! Conceptually, this routine computes <em>u=v+bet*w</em> (if <em>su\f$\neq\f$sv</em>) or <em>v+=bet*w</em> (if <em>su==sv</em>).
 *  <ul>
 *  <li>  <em>v, w</em> are the input vectors of length <em>len</em> stored in <em>sv</em>[], <em>sw</em>[], respectively.
 *  <li>  <em>u</em> is the output vector of length <em>len</em> stored in <em>su</em>[].
 *  </ul>
 *
 *  \param  len       is the length of <em>u, v, w</em>.
 *  \param  sv[]      stores elements of <em>v</em>.
 *  \param  bet       is a constant.
 *  \param  sw[]      stores elements of <em>w</em>.
 *  \param  su[]      stores elements of <em>u</em>.
 *
 *  This routine can be used for computing <em>C=A+bet*B</em> or <em>A+=bet*B</em>,
 *  where <em>A</em> and <em>B</em> are both general or both symmetric matrices.
 *
 *  \remark  If <em>su</em>==<tt>NULL</tt> is given, then the required memory of size <em>len</em> will be allocated for the output vector <em>u</em>;
 *           otherwise it is assumed that sufficient memory has been allocated, with the address pointed to by <em>su</em>.
 *  \remark  This routine invokes the level 1 <b>BLAS</b> routine <b>xAXPY</b> if <tt>USE_BLAS</tt> is defined.
 */
// compute su[i] = sv[i] + bet*sw[i] (if su!=sv) or sw[i] += bet*sv[i] (if su==sv) for i=0,...,len-1
// conceptually, this routine computes u = v+bet*w (if su!=sv) or v += bet*w (if su==sv), where
//    v, w are the input vectors of length len stored in sv[], wv[], respectively, and
//    u is the output vector of length len stored in su[]
// if su==NULL is given, then the required memory of size len will be allocated for the output vector u;
//    otherwise it is assumed that sufficient memory has been allocated, with the address pointed to by su
// this routine can be used for computing C = A+bet*B or A += bet*B, where
//    A and B are both general or both symmetric matrices
// this routine invokes the level 1 BLAS routine xAXPY if USE_BLAS is defined
void elementwise_addition(mkIndex len, const Real *sv, Real bet, const Real *sw, Real *&su);

//! Compute <em>dat_out</em>[<em>perm</em>[<em>i</em>]]=<em>dat_in</em>[<em>i</em>] (if <em>indexFrom1</em>==<b>false</b>) or <em>dat_out</em>[<em>perm</em>[<em>i</em>]<em>-1</em>]<em>=dat_in</em>[<em>i</em>] (if <em>indexFrom1</em>==<b>true</b>) for <em>i=0,...,len-1</em>, where <em>dat_in</em>[] is <em>dat</em>[] as input, and <em>dat_out</em>[] is <em>dat</em>[] as output.
/*!
 *  \param  len       (as input) is the length of the arrays.
 *  \param  perm[]    (as input) is the permutation array.
 *  \param  dat[]     (as input/output) is the the data array.
 *  \param  indexFrom1 (as input) tells whether the indices <em>perm</em>[] start with <em>1</em> or <em>0</em>. \n
 *                    If <em>indexFrom1</em>==<b>true</b>,  the indices <em>perm</em>[] range over <em>1,...,len</em>. \n
 *                    If <em>indexFrom1</em>==<b>false</b>, the indices <em>perm</em>[] range over <em>0,...,len-1</em>.
 *
 *  \return  <em>0</em> on success, or <em>-1</em> if the input <em>perm</em>[] is not a permutation array.
 *
 *  \remark  This is an in-place implementation.
 */
// compute dat_out[perm[i]] = dat_in[i] (if indexFrom1==false) or dat_out[perm[i]-1] = dat_in[i] (if indexFrom1==true) for i=0,...,len-1,
//    where dat_in[] is dat[] as input, and dat_out[] is dat[] as output
// return 0 on success, or -1 if the input perm[] is not a permutation array
// this is an in-place implementation
template<class Index, class Element>
int in_place_permute_elements(mkIndex len, Index *perm, Element *dat, bool indexFrom1=true) {
    // flag all not yet permuted
    for (mkIndex i=0; i<len; i++)
        perm[i] += len;

    // now do the permutation
    Index flag = len + indexFrom1;
    for (Index i=0; i<len; i++) {
        // each inner loop deals with a cycle in the permutation
        if (perm[i] >= flag) {
            Index idx0 = i;
            Element v0 = dat[idx0];
            do {
                perm[idx0] -= len;       // flag entry idx0 for being permuted
                Index idx1 = perm[idx0]-indexFrom1;
                Element v1 = dat[idx1];  // store dat[idx1]
                dat[idx1] = v0;          // move dat[idx0] to dat[idx1]
                v0 = v1;                 // value for the next iteration
                idx0 = idx1;             // index for the next iteration
            } while (perm[idx0] >= flag);
            if (idx0 != i) {
                // the given permutation is not valid!
                return -1;
            }
        }
    }

    // permutation done successfully
    return 0;
}

//! Compute <em>iperm</em>[<em>perm</em>[<em>i</em>]]=<em>i</em> (if <em>indexFrom1</em>==<b>false</b>) or <em>iperm</em>[<em>perm</em>[<em>i</em>]-<em>1</em>]=<em>i+1</em> (if <em>indexFrom1</em>==<b>true</b>) for <em>i=0,...,len-1</em>.
/*!
 *  \param  len       (as input) is the length of the permutation arrays.
 *  \param  perm[]    (as input) is the permutation array.
 *  \param  iperm[]   (as output) is the inverse of the permutation array.
 *  \param  indexFrom1 (as input) tells whether the indices <em>perm</em>[] and <em>iperm</em>[] start with <em>1</em> or <em>0</em>. \n
 *                    If <em>indexFrom1</em>==<b>true</b>,  the indices <em>perm</em>[] and <em>iperm</em>[] range over <em>1,...,len</em>. \n
 *                    If <em>indexFrom1</em>==<b>false</b>, the indices <em>perm</em>[] and <em>iperm</em>[] range over <em>0,...,len-1</em>.
 *
 *  \remark  If <em>iperm</em>==<tt>NULL</tt> is given, then the required memory of size <em>len</em> will be allocated for the inverse permutation;
 *           otherwise it is assumed that sufficient memory has been allocated, with the address pointed to by <em>iperm</em>.
 */
// compute iperm[perm[i]] = i (if indexFrom1==false) or iperm[perm[i]-1] = i+1 (if indexFrom1==true) for i=0,...,len-1
template<class Index>
void inverse_perm(mkIndex len, Index *perm, Index *&iperm, bool indexFrom1=true) {
    if (iperm == NULL) {
        // allocate memory for output
        iperm = new Index[len];
    }

    // shifting if indexFrom1==true
    mkIndex xlen = len + indexFrom1;
    Index *xperm  = perm  - indexFrom1;
    Index *xiperm = iperm - indexFrom1;

    // compute the inverse permutation
    for (Index i=indexFrom1; i<xlen; i++) {
        xiperm[xperm[i]] = i;
    }
}

//! Compute <em>perm_out</em>[<em>perm_in</em>[<em>i</em>]]=<em>i</em> (if <em>indexFrom1</em>==<b>false</b>) or <em>perm_out</em>[<em>perm_in</em>[<em>i</em>]-<em>1</em>]=<em>i+1</em> (if <em>indexFrom1</em>==<b>true</b>) for <em>i=0,...,len-1</em>, where <em>perm_in</em>[] is <em>perm</em>[] as input, and <em>perm_out</em>[] is <em>perm</em>[] as output.
/*!
 *  \param  len       (as input) is the length of the permutation arrays.
 *  \param  perm[]    (as input/output) is the permutation array.
 *  \param  indexFrom1 (as input) tells whether the indices <em>perm</em>[] start with <em>1</em> or <em>0</em>. \n
 *                    If <em>indexFrom1</em>==<b>true</b>,  the indices <em>perm</em>[] range over <em>1,...,len</em>. \n
 *                    If <em>indexFrom1</em>==<b>false</b>, the indices <em>perm</em>[] range over <em>0,...,len-1</em>.
 *
 *  \return  <em>0</em> on success, or <em>-1</em> if the input <em>perm</em>[] is not a permutation array.
 *
 *  \remark  This is an in-place implementation.
 */
// compute perm_out[perm_in[i]] = i (if indexFrom1==false) or perm_out[perm_in[i]-1] = i+1 (if indexFrom1==true) for i=0,...,len-1
//    where perm_in[] is perm[] as input, and perm_out[] is perm[] as output
// this is an in-place implementation
template<class Index>
int in_place_inverse_perm(mkIndex len, Index *perm, bool indexFrom1=true) {
    // flag all as not (yet) permuted
    for (mkIndex i=0; i<len; i++)
        perm[i] += len;

    // now do the permutation
    Index flag = len + indexFrom1;
    for (Index i=0; i<len; i++) {
        // each inner loop deals with a cycle in the permutation
        if (perm[i] >= flag) {
            Index idx0 = i;
            Index idx1 = perm[idx0];
            while (idx1 >= flag) {
                idx1 -= flag;
                Index idx2 = perm[idx1];       // the index, subject to being shifted by flag, of the next entry to update
                perm[idx1] = idx0+indexFrom1;  // update perm[idx1]
                idx0 = idx1;                   // the value, subject to being shifted by indexFrom1, to fill perm[idx2-flag] in the next inner iteration
                idx1 = idx2;                   // the index, subject to being shifted by flag, of the next entry to update
            }
            if (idx1 != i+indexFrom1) {
                // the given permutation is not valid!
                return -1;
            }
        }
    }

    // permutation done successfully
    return 0;
}


#ifdef USE_MEX
// it happens often that an input argument of OCTAVE/MATLAB is a single number
// (e.g., the 2nd input argument of eigs indicates the number of eigenvalues/eigenvectors requested)
// the following routine converts such an input argument to a single real number
// always output "double" even if "Real" is "float", to reduce loss of significant digits due to data-type conversion
double mxArray_to_double(const mxArray *ma, mkIndex offset=0);
#endif

/** @} */  // end of group_misc



////////////////////////////////////////////////////////////////////////////
//    vector routines
////////////////////////////////////////////////////////////////////////////

/** @defgroup group_vector Vector Routines
 *
 *  This module contains routines `without' classes related to vector operations
 *  (declared in matkitfunc.h).
 *
 *  @{
 */

//! Print out vector <em>v</em> to ostream <em>s</em>, where <em>v</em> is stored in <em>sv</em>[] of length <em>len</em>.
// print out vector v to ostream s, where v is stored in sv[] of length len
void print_vector(std::ostream &s, mkIndex len, const Real *sv);

//! Write a vector <em>v</em> to a plain text file <em>fileName</em>[], where <em>v</em> is of length <em>len</em> stored in <em>sv</em>[].
/*! The format of the file:
 *  <ul>
 *  <li>  Each line contains a number.
 *  <li>  The first line has the length of the vector as an integer.
 *  <li>  After that, the (<em>i+1</em>)st has has the <em>i</em>th element of the matrix as a <b>Real</b> number, for <em>i=1,...,len</em>.
 *  </ul>
 *
 *  \param  len       (as input) is the length of <em>v</em>.
 *  \param  sv[]      (as input) stores the elements of <em>v</em>.
 *  \param  fileName[]  (as input) is a character string storing the output matrix-market file name.
 *
 *  \return An integer which indicates the status:
 *  <ul>
 *  <li>  <em>0</em>:  a successful write.
 *  <li>  <em>1</em>:  file open error.
 *  <li>  <em>2</em>:  file close error.
 *  <li>  A negative integer:  a number returned by <b>fprintf</b>() which signals an error.
 *  </ul>
 */
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
int vector_write(mkIndex len, const Real *sv, const char fileName[]);

//! Read a vector from a plain text file <em>fileName</em>[] and store the result in <em>sv</em>[] of length <em>len</em>.
/*! The format of the file:
 *  <ul>
 *  <li>  Each line contains a number.
 *  <li>  The first line has the length of the vector as an integer.
 *  <li>  After that, the (<em>i+1</em>)st has has the <em>i</em>th element of the matrix as a <b>Real</b> number, for <em>i=1,...,len</em>.
 *  </ul>
 *
 *  \param  fileName[]  (as input) is a character string storing the input matrix-market file name.
 *  \param  len       (as output) is the length of <em>v</em>.
 *  \param  sv[]      (as output) stores the elements of <em>v</em>.
 *
 *  \return  An integer which indicates the status:
 *  <ul>
 *  <li>  <em>0</em>:  a successful write.
 *  <li>  <em>1</em>:  file open error.
 *  <li>  <em>2</em>:  file close error.
 *  <li>  A negative integer:  a number returned by <b>fprintf</b>() which signals an error.
 *  </ul>
 */
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
int vector_read(const char fileName[], mkIndex &len, Real *&sv);

//! Compute <em>w=v{op}alp</em> (if <em>sv\f$\neq\f$sw</em>) or <em>v{op}=alp</em> (if <em>sv==sw</em>), where <em>{op}==+,-,*,/</em> if <em>flag==0,1,2,3</em>, respectively.
/*! <ul>
 *  <li>  <em>v</em> is the input  vector of length <em>len</em> stored in <em>sv</em>[].
 *  <li>  <em>w</em> is the output vector of length <em>len</em> stored in <em>sw</em>[].
 *  </ul>
 *
 *  \param  len       (as input) is the length of <em>v, w</em>.
 *  \param  sv[]      (as input) stores the elements of <em>v</em>.
 *  \param  alp       (as input) is a scalar.
 *  \param  sw[]      (as output) stores the elements of <em>w</em>.
 *  \param  flag      (as input) specifies the arithmetic operator. \n
 *
 *  \remark  If <em>sw</em>==<tt>NULL</tt> is given, then the required memory of size <em>len</em> will be allocated for the output vector <em>w</em>;
 *           otherwise, it is assumed that sufficient memory has been allocated, with the address pointed to by <em>sw</em>.
 *  \remark  This routine invokes in_place_scale_elements() if <em>sw==sv</em> and <em>flag==2</em> (i.e. <em>v*=alp</em>).
 *  \remark  This routine can also be used to compute <em>B=A{op}alp</em> (with <em>sv==sw</em>) and <em>A{op}=alp</em> (with <em>sv\f$\neq\f$sw</em>), where
 *           <em>A</em> and <em>B</em> are both general or both symmetric matrices.
 *
 *  \warning  This routine does not catch the division-by-zero exception (i.e. <em>alp</em>==<em>0.0</em> and <em>flag</em>==<em>3</em>).
 *            A division-by-zero will return an <tt>NaN</tt> (not-a-number) for each element.
 */
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
void vector_scalar_operation(mkIndex len, const Real *sv, Real alp, Real *&sw, unsigned flag);

//! Form an (<em>len+|k|</em>)-by-(<em>len+|k|</em>) general matrix <em>A</em> with the <em>k</em>th diagonal formed by vector <em>v</em> and otherwise zero.
/*! <ul>
 *  <li>  <em>v</em> is the input vector of length <em>len</em> stored in <em>sv</em>[].
 *  <li>  <em>A</em> is the output (<em>len+|k|</em>)-by-(<em>len+|k|</em>) general matrix <em>A</em> stored in <em>sa</em>[].
 *  </ul>
 *
 *  To be precise,
 *  if <em>k>=0</em>, then <em>A</em>(<em>j,j+k</em>)=<em>v</em>(<em>j</em>) for <em>j=1,...,len</em> and otherwise zero;
 *  if <em>k<0</em>,  then <em>A</em>(<em>j-k,j</em>)=<em>v</em>(<em>j</em>) for <em>j=1,...,len</em> and otherwise zero.
 *
 *  \param  len       (as input) is the length of <em>v</em>.
 *  \param  sv[]      (as input) stores the elements of <em>v</em>.
 *  \param  sa[]      (as output) stores the elements of <em>A</em>.
 *  \param  k         (as input) is the index of diagonal of <em>A</em>.
 *
 *  \remark  If <em>sa</em>==<tt>NULL</tt> is given, then the required memory of size (<em>len+|k|</em>)*(<em>len+|k|</em>) will be allocated for the output matrix <em>A</em>;
 *           otherwise it is assumed that sufficient memory has been allocated, with the address pointed to by <em>sa</em>.
 *  \remark  This routine resembles <em>A</em>=<b>diag</b>(<em>v,k</em>) in <b>OCTAVE</b>/<b>MATLAB</b>.
 */
// form a (len+|k|)-by-(len+|k|) general matrix A with the k-th diagonal formed by vector v and otherwise zero
// v is the input vector of length len stored in sv[]
// A is the output (len+|k|)-by-(len+|k|) general matrix A stored in sa[], such that
//    if k >= 0, then A(j,j+k)=v(j) for j=1,...,len and otherwise zero, or
//    if k <  0, then A(j-k,j)=v(j) for j=1,...,len and otherwise zero
// if sa==NULL is given, then the required memory of size (len+|k|)*(len+|k|) will be allocated for the output matrix A;
//    otherwise it is assumed that sufficient memory has been allocated, with the address pointed to by sa
// this routine resembles A = diag(v,k) in OCTAVE/MATLAB
void vector_to_diag(mkIndex len, const Real *sv, Real *&sa, mkSignedIndex k=0);

//! Form an (<em>len+|k|</em>)-by-(<em>len+|k|</em>) sparse matrix <em>A</em> with the <em>k</em>th diagonal formed by vector <em>v</em> and otherwise zero.
/*! <ul>
 *  <li>  <em>v</em> is the input vector of length <em>len</em> stored in <em>sv</em>[].
 *  <li>  <em>A</em> is the output (<em>len+|k|</em>)-by-(<em>len+|k|</em>) sparse matrix <em>A</em> in CSC format
 *        stored in <em>store</em>[], <em>rowIdx</em>[], <em>colPtr</em>[].
 *  </ul>
 *
 *  To be precise,
 *  if <em>k>=0</em>, then <em>A</em>(<em>j,j+k</em>)=<em>v</em>(<em>j</em>) for <em>j=1,...,len</em> and otherwise zero;
 *  if <em>k<0</em>,  then <em>A</em>(<em>j-k,j</em>)=<em>v</em>(<em>j</em>) for <em>j=1,...,len</em> and otherwise zero.
 *
 *  \param  len       (as input) is the length of <em>v</em>.
 *  \param  sv[]      (as input) stores the elements of <em>v</em>.
 *  \param  store[]   (as output) stores the elements of <em>A</em> in the order as <em>rowIdx</em>[].
 *  \param  rowIdx[]  (as output) stores the row indices of <em>A</em>.
 *  \param  colPtr[]  (as output) stores the column pointers of <em>A</em>.
 *  \param  k         (as input) is the index of diagonal of <em>A</em>.
 *
 *  \remark  If <em>store</em>==<tt>NULL</tt> (<em>rowIdx</em>==<tt>NULL</tt>, or <em>colPtr</em>==<tt>NULL</tt>) is given, then the required memory will be allocated;
 *           otherwise it is assumed that sufficient memory has been allocated, with the address pointed to by <em>store</em> (<em>rowIdx</em>, or <em>colPtr</em>, respectively).
 *  \remark  This routine resembles <em>A</em>=<b>spdiag</b>(<em>v,k</em>) in <b>OCTAVE</b>, a sparse version of <em>A</em>=<b>diag</b>(<em>v,k</em>) in <b>OCTAVE</b>/<b>MATLAB</b>.
 */
// form a (len+|k|)-by-(len+|k|) sparse matrix A with the k-th diagonal formed by a vector v and otherwise zero
// v is the input vector of length len stored in sv[]
// A is the output (len+|k|)-by-(len+|k|) sparse matrix A stored in store[], rowIdx[], colPtr[], such that
//    if k >= 0, A(j,j+k)=v(j) for j=1,...,len and otherwise zero, or
//    if k <  0, A(j-k,j)=v(j) for j=1,...,len and otherwise zero
// if store==NULL (rowIdx==NULL, or colPtr==NULL) is given, then the required memory will be allocated;
//    otherwise it is assumed that sufficient memory has been allocated, with the address pointed to by store (rowIdx, or colPtr, respectively)
// this routine resembles A = spdiag(v,k) in OCTAVE, a sparse version of diag(v,k) in OCTAVE/MATLAB
void vector_to_spdiag(mkIndex len, const Real *sv, Real *&store, mkIndex *&rowIdx, mkIndex *&colPtr, mkSignedIndex k=0);


#ifdef USE_MEX
// convert *mxArray to a vector
void mxArray_to_vector(const mxArray *mA, mkIndex &len, Real *&sv);
// convert a vector (of length len) to *mxArray (of size len-by-1)
mxArray *vector_to_mxArray(mkIndex len, const Real *sv);
#endif

/** @} */  // end of group_vector



////////////////////////////////////////////////////////////////////////////
//    matrix routines
////////////////////////////////////////////////////////////////////////////

/** @defgroup group_matrix Matrix Routines
 *
 *  This module contains routines `without' classes related to matrix operations
 *  (declared in matkitfunc.h).
 *
 *  @{
 */

//! Print out matrix <em>A</em> to ostream <em>s</em>, where <em>A</em> is of size <em>nr</em>-by-<em>nc</em> stored in <em>sa</em>[].
// print out matrix A to ostream s, where A is of size nr-by-nc stored in sa[]
void print_matrix(std::ostream &s, mkIndex nr, mkIndex nc, const Real *sa);

//! Convert a general matrix <em>A</em> to a sparse matrix <em>B</em>, such that mathematically <em>A==B</em>.
/*! <ul>
 *  <li>  <em>A</em> is the <em>nr</em>-by-<em>nc</em> input general matrix stored in <em>sa</em>[].
 *  <li>  <em>B</em> is the <em>nr</em>-by-<em>nc</em> output sparse matrix stored in <em>sb</em>[], <em>rowIdx</em>[], <em>colPtr</em>[].
 *  </ul>
 *
 *  \param  nr        (as input) is the number of rows of <em>A</em>.
 *  \param  nc        (as input) is the number of columns of <em>A</em>.
 *  \param  sa[]      (as input) stores the elements of <em>A</em>.
 *  \param  sb[]      (as output) stores the elements of <em>B</em> in the order as <em>rowIdx</em>[].
 *  \param  rowIdx[]  (as output) stores the row indices of <em>B</em>.
 *  \param  colPtr[]  (as output) stores the column pointers of <em>B</em>.
 *  \param  patternOnly  (as input) tells whether only the sparsity pattern is considered. \n
 *                    If <em>patternOnly</em>==<b>true</b>, then only the sparsity pattern is considered, and <em>sb</em>[] will not be used.
 *
 *  \remark  If <em>sb</em>==<tt>NULL</tt> (<em>rowIdx</em>==<tt>NULL</tt>, or <em>colPtr</em>==<tt>NULL</tt>) is given, then the required memory will be allocated;
 *           otherwise it is assumed that sufficient memory has been allocated, with the address pointed to by <em>sb</em> (<em>rowIdx</em>, or <em>colPtr</em>, respectively).
 *  \remark  This routine resembles <em>B</em>=<b>sparse</b>(<em>A</em>) in <b>OCTAVE</b>/<b>MATLAB</b>.
 */
// convert a general matrix A to a sparse matrix B, such that mathematically A == B, where
//    A is the nr-by-nc input general matrix stored in sa[], and
//    B is the nr-by-nc output sparse matrix in CSC format stored in sb[], rowIdx[], colPtr[]
// if patternOnly==true (default is false), then only the sparsity pattern is considered, and sb[] will not be used
// if sb==NULL (rowIdx==NULL, or colPtr==NULL) is given, then the required memory will be allocated;
//    otherwise it is assumed that sufficient memory has been allocated, with the address pointed to by sb (rowIdx, or colPtr, respectively)
// this routine resembles B = sparse(A) in OCTAVE/MATLAB
void full_to_sparse(mkIndex nr, mkIndex nc, const Real *sa, Real *&sb, mkIndex *&rowIdx, mkIndex *&colPtr, bool patternOnly=false);

//! Compute a vector <em>v</em> formed by the <em>k</em>th diagonal of a general matrix <em>A</em>.
/*! <ul>
 *  <li>  <em>A</em> is the <em>nr</em>-by-<em>nc</em> input general matrix stored in <em>sa</em>[].
 *  <li>  <em>v</em> is the output vector stored in <em>sv</em>[].
 *  </ul>
 *
 *  To be precise, this routine computes <em>v</em>(<em>i</em>)<em>=A</em>(<em>i,i+k</em>) if <em>k>=0</em>,
 *  or <em>v</em>(<em>i</em>)<em>=A</em>(<em>i-k,i</em>) if <em>k<0</em>.
 *
 *  \param  nr        (as input) is the number of rows of <em>A</em>.
 *  \param  nc        (as input) is the number of columns of <em>A</em>.
 *  \param  sa[]      (as input) stores the elements of <em>A</em>.
 *  \param  sv[]      (as output) stores the elements <em>v</em>.
 *  \param  k         (as input) is the index of diagonal of <em>A</em>.
 *
 *  \return  The length of <em>v</em>.
 *
 *  \remark  If <em>sv</em>==<tt>NULL</tt> is given, then the required memory will be allocated for the output vector <em>v</em>;
 *           otherwise, it is assumed that sufficient memory has been allocated, with the address pointed to by <em>sv</em>.
 *  \remark  This routine resembles <em>v</em>=<b>diag</b>(<em>A</em>,<em>k</em>) in <b>OCTAVE</b>/<b>MATLAB</b>.
 */
// compute a vector v formed by the k-th diagonal of a general matrix A
// to be precise, this routine computes v(i)=A(i,i+k) if k>=0, or v(i)=A(i-k,i) if k<0, where
//    A is the nr-by-nc input general matrix stored in sa[]
//    v is the output vector v stored in sv[] containing the k-th diagonal elements (in order) of A
// the return value is the length of v
// if sv==NULL is given, then the required memory will be allocated for the output vector v;
//    otherwise, it is assumed that sufficient memory has been allocated, with the address pointed to by sv
// this routine resembles v = diag(A,k) in OCTAVE/MATLAB
mkIndex matrix_diag(mkIndex nr, mkIndex nc, const Real *sa, Real *&sv, mkSignedIndex k=0);

//! Transpose a general matrix <em>A</em> and store the result as a general matrix <em>B</em> (not in-place).
/*! To be precise, mathematically <em>B</em>(<em>j,i</em>)=<em>A</em>(<em>i,j</em>) for <em>i=1,...,nr</em> and <em>j=1,...,nc</em>.
 *  <ul>
 *  <li>  <em>A</em> is the <em>nr</em>-by-<em>nc</em> input matrix stored in <em>sa</em>[].
 *  <li>  <em>B</em> is the <em>nc</em>-by-<em>nr</em> output matrix stored in <em>sb</em>[].
 *  </ul>
 *
 *  \param  nr        (as input) is the number of rows of <em>A</em>.
 *  \param  nc        (as input) is the number of columns of <em>A</em>.
 *  \param  sa[]      (as input) stores the elements of <em>A</em>.
 *  \param  sb[]      (as output) stores the elements of <em>B</em>.
 *
 *  \remark  If <em>sb</em>==<tt>NULL</tt> is given, then the required memory of size <em>nr*nc</em> will be allocated for the output general matrix <em>B</em>;
 *           otherwise, it is assumed that sufficient memory has been allocated, with the address pointed to by <em>sb</em>.
 *  \remark  This routine resembles <em>B=A</em>' or equivalently <em>B</em>=<b>transpose</b>(<em>A</em>) in <b>OCTAVE</b>/<b>MATLAB</b>.
 */
// transpose a general matrix A and store the result as a general matrix B (not in-place)
// to be precise, mathematically B(j,i)=A(i,j) for i=1,...,nr and j=1,...,nc, where
//    A is the nr-by-nc input matrix A stored in sa[]
//    B, the transpose of A, is the nc-by-nr output matrix B stored in sb[]
// if sb==NULL is given, then the required memory of size nr*nc will be allocated for the output general matrix B;
//    otherwise, it is assumed that sufficient memory has been allocated, with the address pointed to by sb
// this routine resembles B = A' or equivalently B = transpose(A) in OCTAVE/MATLAB
void matrix_transpose(mkIndex nr, mkIndex nc, const Real *sa, Real *&sb);

//! Compute a general matrix <em>B</em> as the submatrix formed by elements in rows <em>i1,...,i2</em> and columns <em>j1,...,j2</em> of a general matrix <em>A</em>.
/*! To be precise, assume that <em>i1<=i2</em> and <em>j1<=j2</em>.
 *  If <em>indexFrom1</em>==<b>true</b>,  <em>B=A</em>(<em>i1:i2,j1:j2</em>).
 *  If <em>indexFrom1</em>==<b>false</b>, <em>B=A</em>(<em>i1+1:i2+1,j1+1:j2+1</em>).
 *
 *  <ul>
 *  <li>  <em>A</em> is the <em>nr</em>-by-<em>nc</em> input general matrix stored in sa[].
 *  <li>  <em>B</em> is the (<em>|i2-i1|+1</em>)-by-(<em>|j2-j1|+1</em>) output general matrix stored in sb[].
 *  </ul>
 *
 *  \param  nr        (as input) is the number of rows of <em>A</em>.
 *  \param  nc        (as input) is the number of columns of <em>A</em>.
 *  \param  i1,i2     (as input) row indices of <em>A</em>.
 *  \param  j1,j2     (as input) column indices of <em>A</em>.
 *  \param  sa[]      (as input) stores the elements of <em>A</em>.
 *  \param  sb[]      (as output) stores the elements of <em>B</em>.
 *  \param  indexFrom1   (as input) tells whether the input row indices <em>i1,i2</em> and column indices <em>j1,j2</em> start with <em>1</em> or <em>0</em>. \n
 *                    If <em>indexFrom1</em>==<b>true</b>,  the row indices <em>i1,i2</em> and column indices range over <em>1,...,nr</em> and <em>1,...,nc</em>, respectively. \n
 *                    If <em>indexFrom1</em>==<b>false</b>, the row indices <em>i1,i2</em> and column indices range over <em>0,...,nr-1</em> and <em>0,...,nc-1</em>, respectively.
 *
 *  \remark  If <em>sb</em>==<tt>NULL</tt> is given, then the required memory of size (<em>|i2-i1|+1</em>)*(<em>|j2-j1|+1</em>) will be allocated for the output matrix <em>B</em>;
 *           otherwise it is assumed that sufficient memory has been allocated, with the address pointed to by <em>sb</em>.
 *  \remark  This routine allows <em>i1>i2</em> and/or <em>j1>j2</em>.
 *           If <em>i1>i2</em> and <em>j1>j2</em>, then the computed is <em>B=A</em>(<em>i1:-1:i2,j1:-1:j2</em>) if <em>indexFrom1</em>==<b>true</b>,
 *           or <em>B=A</em>(<em>i1+1:-1:i2+1,j1+1:-1:j2+1</em>) if <em>indexFrom1</em>==<b>false</b>.
 */
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
void submatrix_of_matrix(mkIndex nr, mkIndex nc, const Real *sa,
                         mkIndex i1, mkIndex i2, mkIndex j1, mkIndex j2, Real *&sb, bool indexFrom1=true);

//! Permute columns of a general matrix <em>A</em> (not in-place). The result is a general matrix <em>B</em>.
/*! If <em>indexFrom1</em>==<b>true</b>,  <em>B</em>(<em>i,perm</em>[<em>j-1</em>])            = <em>A</em>(<em>i,j</em>) for <em>i=1,...,nr</em> and <em>j=1,...,nc</em>.
 *  If <em>indexFrom1</em>==<b>false</b>, <em>B</em>(<em>i,perm</em>[<em>j-1</em>]+<em>1</em>) = <em>A</em>(<em>i,j</em>) for <em>i=1,...,nr</em> and <em>j=1,...,nc</em>.
 *
 *  <ul>
 *  <li>  <em>A</em> is the <em>nr</em>-by-<em>nc</em> input general matrix stored in sa[].
 *  <li>  <em>B</em> is the <em>nr</em>-by-<em>nc</em> output general matrix stored in sb[].
 *  </ul>
 *
 *  \param  nr        (as input) is the number of rows of <em>A</em>.
 *  \param  nc        (as input) is the number of columns of <em>A</em>.
 *  \param  perm[]    (as input) is the permutation array.
 *  \param  sa[]      (as input) stores the elements of <em>A</em>.
 *  \param  sb[]      (as output) stores the elements of <em>B</em>.
 *  \param  indexFrom1   (as input) tells whether the indices <em>perm</em>[] start with <em>1</em> or <em>0</em>. \n
 *                   If <em>indexFrom1</em>==<b>true</b>,  the indices <em>perm</em>[] range over <em>1,...,nc</em>. \n
 *                   If <em>indexFrom1</em>==<b>false</b>, the indices <em>perm</em>[] range over <em>0,...,nc-1</em>.
 *
 *  \remark  If <em>sb</em>==<tt>NULL</tt> is given, then the required memory of size <em>nr*nc</em> will be allocated for the output matrix <em>B</em>;
 *           otherwise it is assumed that sufficient memory has been allocated, with the address pointed to by <em>sb</em>.
 *  \remark  This routine resembles <em>B</em>(<em>:,perm</em>)=<em>A</em> in <b>OCTAVE</b>/<b>MATLAB</b>.
 */
// permute columns of a general matrix A (not in-place)
// the result is a general matrix B, such that
//    B(i,perm[j-1])   = A(i,j) for i=1,...,nr and j=1,...,nc (if indexFrom1==true), or
//    B(i,perm[j-1]+1) = A(i,j) for i=1,...,nr and j=1,...,nc (if indexFrom1==false), where
//    A and B are the nr-by-nc input and output matrices stored in sa[] and sb[], respectively
// the indices perm[0,...,nc-1] range over 1,...,nc (if indexFrom1==true), or over 0,...,nc-1 (if indexFrom1==false)
// if sb==NULL is given, then the required memory of size nr*nc will be allocated for the output matrix B;
//    otherwise it is assumed that sufficient memory has been allocated, with the address pointed to by sb
// this routine resembles B(:,perm) = A in OCTAVE/MATLAB
void permute_matrix_columns(mkIndex nr, mkIndex nc, const mkIndex *perm, const Real *sa, Real *&sb, bool indexFrom1=true);

//! Permute rows of a general matrix <em>A</em> (not in-place). The result is a general matrix <em>B</em>.
/*! If <em>indexFrom1</em>==<b>true</b>,  <em>B</em>(<em>perm</em>[<em>i-1</em>],<em>j</em>)   = <em>A</em>(<em>i,j</em>) for <em>i=1,...,nr</em> and <em>j=1,...,nc</em>.
 *  If <em>indexFrom1</em>==<b>false</b>, <em>B</em>(<em>perm</em>[<em>i-1</em>]+<em>1,j</em>) = <em>A</em>(<em>i,j</em>) for <em>i=1,...,nr</em> and <em>j=1,...,nc</em>.
 *
 *  <ul>
 *  <li>  <em>A</em> is the <em>nr</em>-by-<em>nc</em> input general matrix stored in sa[].
 *  <li>  <em>B</em> is the <em>nr</em>-by-<em>nc</em> output general matrix stored in sb[].
 *  </ul>
 *
 *  \param  nr        (as input) is the number of rows of <em>A</em>.
 *  \param  nc        (as input) is the number of columns of <em>A</em>.
 *  \param  perm[]    (as input) is the permutation array.
 *  \param  sa[]      (as input) stores the elements of <em>A</em>.
 *  \param  sb[]      (as output) stores the elements of <em>B</em>.
 *  \param  indexFrom1   (as input) tells whether the indices <em>perm</em>[] start with <em>1</em> or <em>0</em>. \n
 *                    If <em>indexFrom1</em>==<b>true</b>,  the indices <em>perm</em>[] range over <em>1,...,nc</em>. \n
 *                    If <em>indexFrom1</em>==<b>false</b>, the indices <em>perm</em>[] range over <em>0,...,nc-1</em>.
 *
 *  \remark  If <em>sb</em>==<tt>NULL</tt> is given, then the required memory of size <em>nr*nc</em> will be allocated for the output matrix <em>B</em>;
 *           otherwise it is assumed that sufficient memory has been allocated, with the address pointed to by <em>sb</em>.
 *  \remark  This routine resembles <em>B</em>(<em>perm,:</em>)=<em>A</em> in <b>OCTAVE</b>/<b>MATLAB</b>.
 */
// permute rows of a general matrix A (not in-place)
// the result is a general matrix B, such that
//    B(perm[i-1],j)   = A(i,j) for i=1,...,nr and j=1,...,nc (if indexFrom1==true), or
//    B(perm[i-1]+1,j) = A(i,j) for i=1,...,nr and j=1,...,nc (if indexFrom1==false), where
//    A and B are the nr-by-nc input and output matrices stored in sa[] and sb[], respectively
// the indices perm[1,...,nr] range over 1,...,nr (if indexFrom1==true), or over 0,...,nr-1 (if indexFrom1==false)
// if sb==NULL is given, then the required memory of size nr*nc will be allocated for the output matrix B;
//    otherwise it is assumed that sufficient memory has been allocated, with the address pointed to by sb
// this routine resembles B(perm,:) = A in OCTAVE/MATLAB
void permute_matrix_rows(mkIndex nr, mkIndex nc, const mkIndex *perm, const Real *sa, Real *&sb, bool indexFrom1=true);


#ifdef USE_MEX
// convert *mxArray to a matrix
void mxArray_to_matrix(const mxArray *mA, mkIndex &nr, mkIndex &nc, Real *&sa);
// convert a matrix to *mxArray
mxArray *matrix_to_mxArray(mkIndex nr, mkIndex nc, const Real *sa);
#endif

/** @} */  // end of group_matrix



////////////////////////////////////////////////////////////////////////////
//    symmetric matrix routines
////////////////////////////////////////////////////////////////////////////

/** @defgroup group_symmatrix Symmetric Matrix Routines
 *
 *  This module contains routines `without' classes related to symmetric matrix operations
 *  (declared in matkitfunc.h).
 *
 *  @{
 */

//! Print out symmetric matrix <em>A</em> to ostream <em>s</em>, where <em>A</em> is of size <em>nrc</em>-by-<em>nrc</em> stored in packed form in <em>sa</em>[].
// print out symmetric matrix A to ostream s, where A is of size nrc-by-nrc stored in packed form in sa[]
void print_symmatrix(std::ostream &s, mkIndex nrc, const Real *sa);

//! Convert a symmetric matrix <em>A</em> to a general matrix <em>B</em>, such that mathematically <em>A==B</em>.
/*! <ul>
 *  <li>  <em>A</em> is the <em>nrc</em>-by-<em>nrc</em> input symmetric matrix stored in packed form in <em>sa</em>[].
 *  <li>  <em>B</em> is the <em>nrc</em>-by-<em>nrc</em> output general matrix stored in <em>sb</em>[].
 *  </ul>
 *
 *  \param  nrc       (as input) is the number of rows/columns of <em>A</em>.
 *  \param  sa[]      (as input) stores the elements of <em>A</em>
 *  \param  sb[]      (as output) stores the elements of <em>B</em>.
 *
 *  \remark  If <em>sb</em>==<tt>NULL</tt> is given, then the required memory of size <em>nrc*nrc</em> will be allocated for the output matrix <em>B</em>;
 *           otherwise it is assumed that sufficient memory has been allocated, with the address pointed to by <em>sb</em>.
 */
// convert a symmetric matrix A to a full matrix B, such that mathematically A == B, where
//    A is the nr-by-nc input symmetric matrix stored in sa[]
//    B is the nr-by-nc output general matrix stored in sb[]
// if sb==NULL is given, then the required memory of size nr*nc will be allocated for the output matrix B;
//    otherwise it is assumed that sufficient memory has been allocated, with the address pointed to by sb
void symmetric_to_general(mkIndex nrc, const Real *sa, Real *&sb);

//! Compute a vector <em>v</em> formed by the <em>k</em>th diagonal of a symmetric matrix <em>A</em>.
/*! <ul>
 *  <li>  <em>A</em> is the <em>nrc</em>-by-<em>nrc</em> input symmetric matrix stored in packed form in <em>sa</em>[].
 *  <li>  <em>v</em> is the output vector stored in <em>sv</em>[].
 *  </ul>
 *
 *  In one word, this routine computes <em>v</em>(<em>i</em>)<em>=A</em>(<em>i+|k|,i</em>) for <em>i=1,...,nrc-|k|</em>.
 *
 *  \param  nrc       (as input) is the number of rows/columns of <em>A</em>.
 *  \param  sa[]      (as input) stores the elements of <em>A</em>.
 *  \param  sv[]      (as output) stores the elements <em>v</em>.
 *  \param  k         (as input) is the index of diagonal of <em>A</em>.
 *
 *  \return  The length of <em>v</em>.
 *
 *  \remark  If <em>sv</em>==<tt>NULL</tt> is given, then the required memory will be allocated for the output vector <em>v</em>;
 *           otherwise, it is assumed that when sufficient memory has been allocated, with the address pointed to by <em>sv</em>.
 *  \remark  This routine resembles <em>v</em>=<b>diag</b>(<em>A</em>,<em>k</em>) in <b>OCTAVE</b>/<b>MATLAB</b>.
 */
// compute a vector v formed by the k-th diagonal of a symmetric matrix A
// in one word, this routine computes v(i)=A(i+|k|,i) for i=1,...,nrc-|k|, where
//    A is the nrc-by-nrc input symmetric matrix stored in packed form in sa[], and
//    v is the output vector stored in sv[]
// the return value is the length of v
// if sv==NULL is given, then the required memory will be allocated for output vector v;
//    otherwise, it is assumed that sufficient memory has been allocated, with the address pointed to by sv
// this routine resembles v = diag(A,k) in OCTAVE/MATLAB
mkIndex symmatrix_diag(mkIndex nrc, const Real *sa, Real *&sv, mkSignedIndex k=0);

//! Compute a vector v formed by the elements <em>i,...,j</em> of row/column <em>k</em> of an <em>nrc</em>-by-<em>nrc</em> symmetric matrix <em>A</em>.
/*! To be precise, assume that <em>i<=j</em>.
 *  If <em>indexFrom1</em>==<b>true</b>,  <em>v=A</em>(<em>i:j,k</em>).
 *  If <em>indexFrom1</em>==<b>false</b>, <em>v=A</em>(<em>i+1:j+1,k+1</em>).
 *
 *  <ul>
 *  <li>  <em>A</em> is the <em>nrc</em>-by-<em>nrc</em> input symmetric matrix stored in packed form in <em>sa</em>[].
 *  <li>  <em>v</em> is a vector of size <em>|j-i|+1</em>.
 *  </ul>
 *
 *  \param  nrc       (as input) is the number of rows/columns of <em>A</em>.
 *  \param  sa[]      (as input) stores the elements of <em>A</em>.
 *  \param  i,j,k     (as input) row indices of <em>A</em>.
 *  \param  sv[]      (as output) stores the elements of <em>v</em>.
 *  \param  incv      (as input) is the index increment of the <em>sv</em>[] elements.
 *  \param  indexFrom1 (as input) tells whether the input indices <em>i,j,k</em> start with <em>1</em> or <em>0</em>. \n
 *                    If <em>indexFrom1</em>==<b>true</b>,  the indices <em>i,j,k</em> range over <em>1,...,nrc</em>. \n
 *                    If <em>indexFrom1</em>==<b>false</b>, the indices <em>i,j,k</em> range over <em>0,...,nr-1</em>.
 *
 *  \remark  If <em>sv</em>==<tt>NULL</tt> is given, then the required memory of size <em>|j-i|+1</em> will be allocated for the output vector <em>v</em>, and
 *           incv will be set as 1 regardless of its input value.
 *           Otherwise it is assumed that all memory accesses via the pointer <em>sv</em> are legitimate.
 *  \remark  This routine allows <em>i>j</em>, in which case
 *           the output is <em>v=A</em>(<em>i:-1:j,k</em>) if <em>indexFrom1</em>==<b>true</b>, or <em>B=A</em>(<em>i+1:-1:j+1,k+1</em>) if <em>indexFrom1</em>==<b>false</b>.
 */
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
void elements_of_symmatrix(mkIndex nrc, const Real *sa, mkIndex k, mkIndex i, mkIndex j, Real *&sv, mkSignedIndex incv=1, bool indexFrom1=true);

//! Compute a general matrix <em>B</em> as the submatrix formed by elements in rows <em>i1,...,i2</em> and columns <em>j1,...,j2</em> of a symmetric matrix <em>A</em>.
/*! To be precise, assume that <em>i1<=i2</em> and <em>j1<=j2</em>.
 *  If <em>indexFrom1</em>==<b>true</b>,  <em>B=A</em>(<em>i1:i2,j1:j2</em>).
 *  If <em>indexFrom1</em>==<b>false</b>, <em>B=A</em>(<em>i1+1:i2+1,j1+1:j2+1</em>).
 *
 *  <ul>
 *  <li>  <em>A</em> is the <em>nrc</em>-by-<em>nrc</em> input symmetric matrix stored in packed form in <em>sa</em>[].
 *  <li>  <em>B</em> is the (<em>|i2-i1|+1</em>)-by-(<em>|j2-j1|+1</em>) output general matrix stored in <em>sb</em>[].
 *  </ul>
 *
 *  \param  nrc       (as input) is the number of rows/columns of <em>A</em>.
 *  \param  sa[]      (as input) stores the elements of <em>A</em>.
 *  \param  i1,i2     (as input) row indices of <em>A</em>.
 *  \param  j1,j2     (as input) column indices of <em>A</em>.
 *  \param  sb[]      (as output) stores the elements of <em>B</em>.
 *  \param  indexFrom1 (as input) tells whether the input row indices <em>i1,i2</em> and column indices <em>j1,j2</em> start with <em>1</em> or <em>0</em>. \n
 *                    If <em>indexFrom1</em>==<b>true</b>,  the row indices <em>i1,i2</em> and column indices <em>j1,j2</em> range over <em>1,...,nr</em> and <em>1,...,nc</em>, respectively. \n
 *                    If <em>indexFrom1</em>==<b>false</b>, the row indices <em>i1,i2</em> and column indices <em>j1,j2</em> range over <em>0,...,nr-1</em> and <em>0,...,nc-1</em>, respectively.
 *
 *  \remark  If <em>sb</em>==<tt>NULL</tt> is given, then the required memory of size (<em>|i2-i1|+1</em>)*(<em>|j2-j1|+1</em>) will be allocated for the output matrix <em>B</em>;
 *           otherwise it is assumed that sufficient memory has been allocated, with the address pointed to by <em>sb</em>.
 *  \remark  This routine allows <em>i1>i2</em> and/or <em>j1>j2</em>.
 *           If <em>i1>i2</em> and <em>j1>j2</em>, then the computed is <em>B=A</em>(<em>i1:-1:i2,j1:-1:j2</em>) if <em>indexFrom1</em>==<b>true</b>,
 *           or <em>B=A</em>(<em>i1+1:-1:i2+1,j1+1:-1:j2+1</em>) if <em>indexFrom1</em>==<b>false</b>.
 */
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
void submatrix_of_symmatrix(mkIndex nrc, const Real *sa, mkIndex i1, mkIndex i2, mkIndex j1, mkIndex j2, Real *&sb, bool indexFrom1=true);

//! Compute a symmetric matrix <em>B</em> as the submatrix formed by elements in rows/columns <em>k1,...,k2</em> of a symmetric matrix <em>A</em>.
/*! To be precise, assume that <em>k1<=k2</em>.
 *  If <em>indexFrom1</em>==<b>true</b>,  <em>B=A</em>(<em>k1:k2,k1:k2</em>).
 *  If <em>indexFrom1</em>==<b>false</b>, <em>B=A</em>(<em>k1+1:k2+1,k1+1:k2+1</em>).
 *
 *  <ul>
 *  <li>  <em>A</em> is the <em>nrc</em>-by-<em>nrc</em> input symmetric matrix stored in packed form in <em>sa</em>[].
 *  <li>  <em>B</em> is the (<em>|k2-k1|+1</em>)-by-(<em>|k2-k1|+1</em>) output symmetric matrix stored in packed form in <em>sb</em>[].
 *  </ul>
 *
 *  \param  nrc       (as input) is the number of rows/columns of <em>A</em>.
 *  \param  sa[]      (as input) stores the elements of <em>A</em>.
 *  \param  k1,k2     (as input) row/column indices of <em>A</em>.
 *  \param  sb[]      (as output) stores the elements of <em>B</em>.
 *  \param  indexFrom1 (as input) tells whether the input row/column indices <em>k1,k2</em> start with <em>1</em> or <em>0</em>. \n
 *                    If <em>indexFrom1</em>==<b>true</b>,  the indices <em>k1,k2</em> range over <em>1,...,nr</em> and <em>1,...,nc</em>, respectively. \n
 *                    If <em>indexFrom1</em>==<b>false</b>, the indices <em>k1,k2</em> range over <em>0,...,nr-1</em> and <em>0,...,nc-1</em>, respectively.
 *
 *  \remark  If <em>sb</em>==<tt>NULL</tt> is given, then the required memory of size (<em>|i2-i1|+1</em>)*(<em>|j2-j1|+1</em>) will be allocated for the output matrix <em>B</em>;
 *           otherwise it is assumed that sufficient memory has been allocated, with the address pointed to by <em>store2</em> (<em>rowIdx2</em>, <em>colPtr2</em>, respectively).
 *  \remark  This routine allows <em>k1>k2</em>, in which case
 *           the computed is <em>B=A</em>(<em>k1:-1:k2,k1:-1:k2</em>) if <em>indexFrom1</em>==<b>true</b>,
 *           or <em>B=A</em>(<em>k1+1:-1:k2+1,k1+1:-1:k2+1</em>) if <em>indexFrom1</em>==<b>false</b>.
 */
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
void submatrix_of_symmatrix(mkIndex nrc, const Real *sa, mkIndex k1, mkIndex k2, Real *&sb, bool indexFrom1=true);

//! Symmetrically permute rows and columns of a symmetric matrix A (not in-place). The result is a symmetric matrix B.
/*! If <em>indexFrom1</em>==<b>true</b>,  <em>B</em>(<em>perm</em>[<em>i-1</em>],<em>perm</em>[<em>j-1</em>])            = <em>A</em>(<em>i,j</em>) for <em>i,j=1,...,nrc</em>.
 *  If <em>indexFrom1</em>==<b>false</b>, <em>B</em>(<em>perm</em>[<em>i-1</em>]+<em>1</em>,<em>perm</em>[<em>j-1</em>]+<em>1</em>) = <em>A</em>(<em>i,j</em>) for <em>i,j=1,...,nrc</em>.
 *
 *  <ul>
 *  <li>  <em>A</em> is the <em>nrc</em>-by-<em>nrc</em> input  symmetric matrix stored in packed form in <em>sa</em>[].
 *  <li>  <em>B</em> is the <em>nrc</em>-by-<em>nrc</em> output symmetric matrix stored in packed form in <em>sb</em>[].
 *  </ul>
 *
 *  \param  nrc       (as input) is the number of rows/columns of <em>A</em>.
 *  \param  perm[]    (as input) is the permutation array.
 *  \param  sa[]      (as input) stores the elements of <em>A</em> in packed form.
 *  \param  sb[]      (as input) stores the elements of <em>B</em> in packed form.
 *  \param  indexFrom1 (as input) tells whether the indices <em>perm</em>[] start with <em>1</em> or <em>0</em>. \n
 *                    If <em>indexFrom1</em>==<b>true</b>,  the indices <em>perm</em>[] range over <em>1,...,nrc</em>. \n
 *                    If <em>indexFrom1</em>==<b>false</b>, the indices <em>perm</em>[] range over <em>0,...,nrc-1</em>.
 *
 *  \remark  If <em>sb</em>==<tt>NULL</tt> is given, then the required memory of size <em>nrc</em>*(<em>nrc+1</em>)/<em>2</em> will be allocated for the output symmetric matrix <em>B</em>;
 *           otherwise it is assumed that sufficient memory has been allocated, with the address pointed to by <em>sb</em>.
 *  \remark  This routine resembles <em>B</em>(<em>perm,perm</em>)=<em>A</em> in <b>OCTAVE</b>/<b>MATLAB</b>.
 */
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
void permute_symmatrix_rows_and_columns(mkIndex nrc, const mkIndex *perm, const Real *sa, Real *&sb, bool indexFrom1=true);

/** @} */  // end of group_symmatrix



////////////////////////////////////////////////////////////////////////////
//    sparse matrix routines
////////////////////////////////////////////////////////////////////////////

/** @defgroup group_spmatrix Sparse Matrix Routines
 *
 *  This module contains routines `without' classes related to sparse matrix operations
 *  (declared in matkitfunc.h).
 *
 *  @{
 */

//! Print out sparse matrix <em>A</em> to ostream <em>s</em>, where <em>A</em> is of size <em>nr</em>-by-<em>nc</em> in CSC format stored in <em>sa</em>[], <em>rowIdx</em>[], <em>colPtr</em>[].
/*! \remark  If indexFrom1==true (or indexFrom1==false), the printed row/column indices start from 1 (or 0, respectively). 
 *  \remark  If <em>patternOnly</em>==<b>true</b>, then only the sparsity pattern is to be printed, and <em>sa</em>[] will not be used.
 */
// print out sparse matrix A to ostream s, where A is of size nr-by-nc in CSC format stored in sa[], rowIdx[], colPtr[]
// if indexFrom1==true (or indexFrom1==false), the printed row/column indices start from 1 (or 0, respectively)
// if patternOnly==true, then only the sparsity pattern is to be printed, and sa[] will not be used
void print_spmatrix(std::ostream &s, mkIndex nr, mkIndex nc, const Real *sa, const mkIndex *rowIdx, const mkIndex *colPtr, bool indexFrom1=true, bool patternOnly=false);

//! Write a sparse matrix <em>A</em> to a matrix-market file <em>fileName</em>[], where <em>A</em> is in CSC format stored in <em>sa</em>[], <em>rowIdx</em>[], <em>colPtr</em>[].
/*! \param  nrows     (as input) is the number of rows of <em>A</em>.
 *  \param  ncols     (as input) is the number of columns of <em>A</em>.
 *  \param  sa[]      (as input) stores the nonzero elements of <em>A</em> in the order as <em>rowIdx</em>[].
 *  \param  rowIdx[]  (as input) stores the row indices of <em>A</em>.
 *  \param  colPtr[]  (as input) stores the column pointers of <em>A</em>.
 *  \param  fileName[]  (as input) is a character string storing the output matrix-market file name.
 *  \param  patternOnly  (as input) tells whether only the sparsity pattern is considered. \n
 *                    If <em>patternOnly</em>==<b>true</b>, then only the sparsity pattern is to be written to the file <em>fileName</em>[], and <em>sa</em>[] will not be used.
 *  \param  comment[]  (as input) is a character string as a comment. \n
 *                    If \f$\neq\f$<tt>NULL</tt>, then the comment will be printed in the output matrix-market file.
 *
 *  \return  An integer which indicates the status:
 *  <ul>
 *  <li>  <em>0</em>:  a successful write.
 *  <li>  <em>1</em>:  file open error.
 *  <li>  <em>2</em>:  file close error.
 *  <li>  A negative integer:  a number returned by <b>fprintf</b>() which signals an error.
 *  </ul>
 */
// write a sparse matrix A to a matrix-market file fileName[],
//    where A is in CSC format stored in sa[], rowIdx[], colPtr[]
// nrows, ncols are the number of rows, the number of columns, respectively
// if patternOnly==true, then only the sparsity pattern is to be written to the file fileName[], and
//    sa[] will not be used
// comment[] is a character string which will be printed in the output file fileName[] (if comment!=NULL)
// the return value:
//    0: a successful read
//    1: file open error
//    2: file close error
//    a negative integer: a number returned by fprintf() which signals an error
int matrix_market_spmatrix_write(mkIndex nrows, mkIndex ncols, const Real *sa, const mkIndex *rowIdx, const mkIndex *colPtr,
                                 const char fileName[], bool patternOnly=false, const char comment[]=NULL);

//! Read a sparse matrix <em>A</em> from a matrix-market file <em>fileName</em>[] and store the result in the coordinate format in <em>sa</em>[], <em>ridx</em>[], <em>cidx</em>[].
/*! \param  fileName[]  (as input) is a character string storing the input matrix-market file name.
 *  \param  nrows     (as output) is the number of rows of <em>A</em>.
 *  \param  ncols     (as output) is the number of columns of <em>A</em>.
 *  \param  nnz       (as output) is the number of nonzero elements of <em>A</em>.
 *  \param  sa[]      (as output) stores the nonzero elements of <em>A</em> in the order as <em>ridx</em>[], <em>cidx</em>[].
 *  \param  ridx[]    (as output) stores the row indices of <em>A</em>.
 *  \param  cidx[]    (as output) stores the column indices of <em>A</em>..
 *  \param  symm      (as output) gives the symmetry information: 'g' for general, 's' for symmetric, and 'k' for skew-symmetric. \n
 *                    If <em>A</em> is symmetric or skew-symmetric, then each pair of nonzero off-diagonal elements <em>A</em>(<em>i,j</em>) and <em>A</em>(<em>j,i</em>) is recorded once.
 *
 *  \return  An integer which indicates the status:
 *  <ul>
 *  <li>  <em>0</em>:  a successful read.
 *  <li>  <em>1</em>:  the input file is a dense matrix.
 *  <li>  <em>2</em>:  the input file is a complex (sparse) matrix.
 *  <li>  <em>3</em>:  the input file contains only sparsity information.
 *  <li>  <em>-1</em>: file open error.
 *  <li>  <em>-2</em>: file is empty.
 *  <li>  <em>-3</em>: the first line (the header) of the file is empty.
 *  <li>  <em>-4</em>: invalid header.
 *  <li>  <em>-5</em>: invalid number of rows, number of columns, or number of nonzero elements.
 *  <li>  <em>-6</em>: the matrix data is incomplete or invalid.
 *  </ul>
 *
 *  \remark  If <em>sa</em>==<tt>NULL</tt> (<em>ridx</em>==<tt>NULL</tt>, or <em>cidx</em>==<tt>NULL</tt>) is given, then the required memory will be allocated;
 *           otherwise it is assumed that sufficient memory has been allocated, with the address pointed to by <em>sa</em> (<em>ridx</em>, <em>cidx</em>, respectively)
 *  \remark  It also takes matrix files in the straight coordinate format (i.e., without the matrix-market header), in which case <em>symm</em>='<tt>g</tt>' will be set.
 */
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
int matrix_market_spmatrix_read(const char fileName[], mkIndex &nrows, mkIndex &ncols, mkIndex &nnz,
                                Real *&sa, mkIndex *&ridx, mkIndex *&cidx, char &symm);

//! Convert a sparse matrix <em>A</em> in coordinate format to the compressed sparse column (CSC) format.
/*! <ul>
 *  <li>  The input coordinate format is specified by <em>sa</em>[], <em>ridx</em>[], <em>cidx</em>[].
 *  <li>  The output CSC format is specified by <em>store</em>[], <em>rowIdx</em>[], <em>colPtr</em>[].
 *  </ul>
 *
 *  \param  nc        (as input) is the number of columns of <em>A</em>.
 *  \param  nnz       (as input) is the number of nonzero elements of <em>A</em>.
 *  \param  sa[]      (as input) stores the elements of <em>A</em> in the order as <em>ridx</em>[], <em>cidx</em>[].
 *  \param  ridx[]    (as input) stores the row indices of <em>A</em>.
 *  \param  cidx[]    (as input) stores the column indices of <em>A</em>.
 *  \param  store[]   (as output) stores the elements of <em>A</em> in the order as <em>rowIdx</em>[].
 *  \param  rowIdx[]  (as output) stores the row indices of <em>A</em>.
 *  \param  colPtr[]  (as output) stores the column pointers of <em>A</em>.
 *  \param  symm      (as input) gives the symmetry information: 'g' for general, 's' for symmetric, and 'k' for skew-symmetric. \n
 *                    If <em>A</em> is symmetric or skew-symmetric, then each pair of nonzero off-diagonal elements <em>A</em>(<em>i,j</em>) and <em>A</em>(<em>j,i</em>) is recorded once in the coordinate format (as input);
 *                    however the elements will be duplicated in the output CSC format (as output) for being symmetric or skew-symmetric.
 *  \param  indexFrom1  (as input) tells whether the input row indices <em>ridx</em>[] and column indices <em>cidx</em>[] start with <em>1</em> or <em>0</em>. \n
 *                    If <em>indexFrom1</em>==<b>true</b>,  the row indices <em>ridx</em>[] and the column indices <em>cidx</em>[] range over <em>1,...,nr</em> and <em>1,...,nc</em>, respectively. \n
 *                    If <em>indexFrom1</em>==<b>false</b>, the row indices <em>ridx</em>[] and the column indices <em>cidx</em>[] range over <em>0,...,nr-1</em> and <em>0,...,nc-1</em>, respectively. \n
 *                    The output row indices <em>rowIdx</em>[] always range over <em>0,...,nr-1</em>.
 *  \param  patternOnly  (as input) tells whether only the sparsity pattern is considered. \n
 *                    If <em>patternOnly</em>==<b>true</b>, then only the sparsity pattern is considered, and <em>sa</em>[] and <em>store</em>[] will not be used.
 *
 *  \return  <em>0</em> on success, or <em>-1</em> if <em>symm</em> is not valid.
 *
 *  \remark  The number of rows of <em>A</em>, <em>nr</em>, is not required in the computation, so it is not passed.
 *  \remark  If <em>store</em>==<tt>NULL</tt> (<em>rowIdx</em>==<tt>NULL</tt>, or <em>colPtr</em>==<tt>NULL</tt>) is given, then the required memory will be allocated;
 *           otherwise it is assumed that sufficient memory has been allocated, with the address pointed to by <em>store</em> (<em>rowIdx</em>, <em>colPtr</em>, respectively).
 *  \remark  This routine can also convert a sparse matrix in coordinate format to the compressed sparse row (CSR) format,
 *           by swapping <em>ridx</em>[] and <em>cidx</em>[] and passing <em>nc</em> as the number of rows.
 */
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
int coordinates_to_csc(mkIndex nc, mkIndex nnz, const Real *sa, const mkIndex *ridx, const mkIndex *cidx,
                       Real *&store, mkIndex *&rowIdx, mkIndex *&colPtr,
                       char symm='g', bool indexFrom1=true, bool patternOnly=false);

//! Sort elements in each column of a sparse matrix <em>A</em> in CSC format (in-place implementation).
/*! <ul>
 *  <li>  <em>A</em> is the <em>nr</em>-by-<em>nc</em> input sparse matrix in CSC format stored in <em>store</em>[], <em>rowIdx</em>[], <em>colPtr</em>[],
 *        but the elements in each column can be in any order (i.e. unsorted).
 *  <li>  The output is the same matrix in CSC format, using the input storage <em>store</em>[], <em>rowIdx</em>[], <em>colPtr</em>[],
 *        with the elements in each column sorted with respect to the row indices.
 *  </ul>
 *
 *  \param  nr        (as input) is the number of rows of <em>A</em>.
 *  \param  nc        (as input) is the number of columns of <em>A</em>.
 *  \param  store[]   (as input/output) stores the elements of <em>A</em> in the order as <em>rowIdx</em>[].
 *  \param  rowIdx[]  (as input/output) stores the row indices of <em>A</em>.
 *  \param  colPtr[]  (as input/output) stores the column pointers of <em>A</em>.
 *  \param  patternOnly  (as input) tells whether only the sparsity pattern is considered. \n
 *                    If <em>patternOnly</em>==<b>true</b>, then only the sparsity pattern is considered, and <em>store</em>[] will not be used.
 *  \param  work[]    (as input) is a work space, which should be of size at least <em>nnz</em>+<b>max</b>(<em>nnz,nr+1</em>).
 *
 *  \return  <em>0</em> on success, or <em>-1</em> if the computed permutation is not valid.
 *
 *  \remark  If <em>work</em>==<tt>NULL</tt> is given, then work space allowing <em>nnz</em>+<b>max</b>(<em>nnz,nr+1</em>) <b>mkIndex</b>s is allocated for sorting, and
 *           the space will be freed when sorting is done.
 *           Otherwise, it is assumed that sufficient memory has been allocated with the address pointed to by <em>work</em>, and it will not be freed in this routine.
 *  \remark  This routine can be used to sort elements in each row of a sparse matrix in CSR format.
 */
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
int sort_csc_elements_in_each_column(mkIndex nr, mkIndex nc, Real *store, mkIndex *rowIdx, mkIndex *colPtr, bool patternOnly=false, mkIndex *work=NULL);

//! Convert a sparse matrix <em>A</em> to a full matrix <em>B</em>, such that mathematically <em>A==B</em>.
/*! <ul>
 *  <li>  <em>A</em> is the <em>nr</em>-by-<em>nc</em> input sparse matrix stored in <em>sa</em>[], <em>rowIdx</em>[], <em>colPtr</em>[].
 *  <li>  <em>B</em> is the <em>nr</em>-by-<em>nc</em> output general matrix stored in <em>sb</em>[].
 *  </ul>
 *
 *  \param  nr        (as input) is the number of rows of <em>A</em>.
 *  \param  nc        (as input) is the number of columns of <em>A</em>.
 *  \param  sa[]      (as input) stores the elements of <em>A</em> in the order as <em>rowIdx</em>[].
 *  \param  rowIdx[]  (as input) stores the row indices of <em>A</em>.
 *  \param  colPtr[]  (as input) stores the column pointers of <em>A</em>.
 *  \param  sb[]      (as output) stores the elements of <em>B</em>.
 *
 *  \remark  If <em>sb</em>==<tt>NULL</tt> is given, then the required memory of size <em>nr*nc</em> will be allocated for the output matrix <em>B</em>;
 *           otherwise it is assumed that sufficient memory has been allocated, with the address pointed to by <em>sb</em>.
 *  \remark  This routine resembles <em>B</em>=<b>full</b>(<em>A</em>) in <b>OCTAVE</b>/<b>MATLAB</b>.
 */
// convert a sparse matrix A to a full matrix B, such that mathematically A == B, where
//    A is the nr-by-nc input sparse matrix stored in sa[], rowIdx[], colPtr[]
//    B is the nr-by-nc output general matrix stored in sb[]
// if sb==NULL is given, then the required memory of size nr*nc will be allocated for the output matrix B;
//    otherwise it is assumed that sufficient memory has been allocated, with the address pointed to by sb
// this routine resembles B = full(A) in OCTAVE/MATLAB
void sparse_to_full(mkIndex nr, mkIndex nc, const Real *sa, const mkIndex *rowIdx, const mkIndex *colPtr, Real *&sb);

//! Compute a vector <em>v</em> formed by the <em>k</em>th diagonal of a sparse matrix <em>A</em>.
/*! <ul>
 *  <li>  <em>A</em> is the <em>nr</em>-by-<em>nc</em> input sparse matrix in CSC format stored in <em>store</em>[], <em>rowIdx</em>[], <em>colPtr</em>[].
 *  <li>  <em>v</em> is the output vector stored in <em>sv</em>[].
 *  </ul>
 *
 *  To be precise, this routine computes <em>v</em>(<em>i</em>)<em>=A</em>(<em>i,i+k</em>) if <em>k>=0</em>,
 *  or <em>v</em>(<em>i</em>)<em>=A</em>(<em>i-k,i</em>) if <em>k<0</em>.
 *
 *  \param  nr        (as input) is the number of rows of <em>A</em>.
 *  \param  nc        (as input) is the number of columns of <em>A</em>.
 *  \param  store[]   (as input) stores the elements of <em>A</em> in the order as <em>rowIdx</em>[].
 *  \param  rowIdx[]  (as input) stores the row indices of <em>A</em>.
 *  \param  colPtr[]  (as input) stores the column pointers of <em>A</em>.
 *  \param  sv[]      (as output) stores the elements <em>v</em>.
 *  \param  k         (as input) is the index of diagonal of <em>A</em>.
 *
 *  \return  The length of <em>v</em>.
 *
 *  \remark  If <em>sv</em>==<tt>NULL</tt> is given, then the required memory will be allocated for the output vector v;
 *           otherwise, it is assumed that sufficient memory has been allocated, with the address pointed to by <em>sv</em>.
 *  \remark  This routine can also be used to extract the k-th diagonal of a sparse matrix in CSR format, by changing the sign of <em>k</em>.
 *  \remark  This routine resembles <em>v</em>=<b>diag</b>(<em>A</em>,<em>k</em>) in <b>OCTAVE</b>/<b>MATLAB</b>.
 */
// compute a vector v formed by the k-th diagonal of a sparse matrix A
// to be precise, this routine computes v(i)=A(i,i+k) if k>=0, or v(i)=A(i-k,i) if k<0, where
//    A is the nr-by-nc input sparse matrix in CSC format stored in store[], rowIdx[], colPtr[]
//    v is the output vector v stored in sv[] containing the k-th diagonal elements (in order) of A
// the return value is the length of v
// if sv==NULL is given, then the required memory will be allocated for the output vector v;
//    otherwise, it is assumed that sufficient memory has been allocated, with the address pointed to by sv
// this routine can also be used to find the k-th diagonal of a sparse matrix in CSR format, by changing the sign of k
// this routine resembles v = diag(A,k) in OCTAVE/MATLAB
mkIndex spmatrix_diag(mkIndex nr, mkIndex nc, const Real *store, const mkIndex *rowIdx, const mkIndex *colPtr, Real *&sv, mkSignedIndex k=0);

//! Transpose a sparse matrix <em>A</em> and store the result as a sparse matrix <em>B</em> (not in-place).
/*! To be precise, mathematically <em>B</em>(<em>j,i</em>)=<em>A</em>(<em>i,j</em>) for <em>i=1,...,nr</em> and <em>j=1,...,nc</em>.
 *  <ul>
 *  <li>  <em>A</em> is the <em>nr</em>-by-<em>nc</em> input sparse matrix in CSC format stored in <em>store</em>[], <em>rowIdx</em>[], <em>colPtr</em>[].
 *  <li>  <em>B</em> is the <em>nc</em>-by-<em>nr</em> output sparse matrix in CSC format stored in <em>store2</em>[], <em>rowIdx2</em>[], <em>colPtr2</em>[].
 *  </ul>
 *
 *  \param  nr        (as input) is the number of rows of <em>A</em>.
 *  \param  nc        (as input) is the number of columns of <em>A</em>.
 *  \param  store[]   (as input) stores the elements of <em>A</em> in the order as <em>rowIdx</em>[].
 *  \param  rowIdx[]  (as input) stores the row indices of <em>A</em>.
 *  \param  colPtr[]  (as input) stores the column pointers of <em>A</em>.
 *  \param  store2[]  (as output) stores the elements of <em>B</em> in the order as <em>rowIdx2</em>[].
 *  \param  rowIdx2[] (as output) stores the row indices of <em>B</em>.
 *  \param  colPtr2[] (as output) stores the column pointers of <em>B</em>.
 *  \param  patternOnly  (as input) tells whether only the sparsity pattern is considered. \n
 *                    If <em>patternOnly</em>==<b>true</b>, then only the sparsity pattern is considered, and <em>store</em>[] and <em>store2</em>[] will not be used.
 *
 *  \remark  If <em>store2</em>==<tt>NULL</tt> (<em>rowIdx2</em>==<tt>NULL</tt>, or <em>colPtr2</em>==<tt>NULL</tt>) is given, then the required memory will be allocated;
 *           otherwise it is assumed that sufficient memory has been allocated, with the address pointed to by <em>store2</em> (<em>rowIdx2</em>, or <em>colPtr2</em>, respectively).
 *  \remark  The output sparse matrix <em>B</em> has elements in each column sorted with respect to row indices,
 *           no matter whether the input matrix has this property or not.
 *  \remark  This routine resembles <em>B=A</em>' or equivalently <em>B</em>=<b>transpose</b>(<em>A</em>) in <b>OCTAVE</b>/<b>MATLAB</b>.
 *  \remark  This routine can also be used for transposing a sparse matrix in CSR format.
 */
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
                        Real *&store2, mkIndex *&rowIdx2, mkIndex *&colPtr2, bool patternOnly=false);

//! Compute a sparse matrix <em>B</em> as the submatrix formed by elements in rows <em>i1,...,i2</em> and columns <em>j1,...,j2</em> of a sparse matrix <em>A</em>.
/*! To be precise, assume that <em>i1<=i2</em> and <em>j1<=j2</em>.
 *  If <em>indexFrom1</em>==<b>true</b>,  <em>B=A</em>(<em>i1:i2,j1:j2</em>).
 *  If <em>indexFrom1</em>==<b>false</b>, <em>B=A</em>(<em>i1+1:i2+1,j1+1:j2+1</em>).
 *
 *  <ul>
 *  <li>  <em>A</em> is the <em>nr</em>-by-<em>nc</em> input sparse matrix in CSC format stored in <em>store</em>[], <em>rowIdx</em>[], <em>colPtr</em>[].
 *  <li>  <em>B</em> is the (<em>|i2-i1|+1</em>)-by-(<em>|j2-j1|+1</em>) output sparse matrix in CSC format stored in <em>store2</em>[], <em>rowIdx2</em>[], <em>colPtr2</em>[].
 *  </ul>
 *
 *  \param  nc        (as input) is the number of columns of <em>A</em>.
 *  \param  store[]   (as input) stores the elements of <em>A</em> in the order as <em>rowIdx</em>[].
 *  \param  rowIdx[]  (as input) stores the row indices of <em>A</em>.
 *  \param  colPtr[]  (as input) stores the column pointers of <em>A</em>.
 *  \param  i1,i2     (as input) row indices of <em>A</em>.
 *  \param  j1,j2     (as input) column indices of <em>A</em>.
 *  \param  store2[]  (as output) stores the elements of <em>B</em> in the order as <em>rowIdx2</em>[].
 *  \param  rowIdx2[] (as output) stores the row indices of <em>B</em>.
 *  \param  colPtr2[] (as output) stores the column pointers of <em>B</em>.
 *  \param  indexFrom1   (as input) tells whether the input row indices <em>i1,i2</em> and column indices <em>j1,j2</em> start with <em>1</em> or <em>0</em>. \n
 *                    If <em>indexFrom1</em>==<b>true</b>,  the row indices <em>i1,i2</em> and column indices range over <em>1,...,nr</em> and <em>1,...,nc</em>, respectively. \n
 *                    If <em>indexFrom1</em>==<b>false</b>, the row indices <em>i1,i2</em> and column indices range over <em>0,...,nr-1</em> and <em>0,...,nc-1</em>, respectively.
 *  \param  patternOnly  (as input) tells whether only the sparsity pattern is considered. \n
 *                    If <em>patternOnly</em>==<b>true</b>, then only the sparsity pattern is considered, and <em>store</em>[] and <em>store2</em>[] will not be used.
 *
 *  \remark  The number of rows of <em>A</em>, <em>nr</em>, is not required in the computation, so it is not passed.
 *  \remark  If <em>store2</em>==<tt>NULL</tt> (<em>rowIdx2</em>==<tt>NULL</tt>, or <em>colPtr2</em>==<tt>NULL</tt>) is given, then the required memory will be allocated;
 *           otherwise it is assumed that sufficient memory has been allocated, with the address pointed to by <em>store2</em> (<em>rowIdx2</em>, <em>colPtr2</em>, respectively).
 *  \remark  This routine allows <em>i1>i2</em> and/or <em>j1>j2</em>.
 *           If <em>i1>i2</em> and <em>j1>j2</em>, then the computed is <em>B=A</em>(<em>i1:-1:i2,j1:-1:j2</em>) if <em>indexFrom1</em>==<b>true</b>,
 *           or <em>B=A</em>(<em>i1+1:-1:i2+1,j1+1:-1:j2+1</em>) if <em>indexFrom1</em>==<b>false</b>.
 *  \remark  This routine does not require the input sparse matrix with elements in each column sorted by the row indices.
 *           However, the output matrix <em>B</em> is guaranteed to have this property, only if the input matrix <em>A</em> has this property.
 */
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
                           Real *&store2, mkIndex *&rowIdx2, mkIndex *&colPtr2, bool indexFrom1=true, bool patternOnly=false);

//! Permute columns of a sparse matrix <em>A</em> in CSC format (not in-place). The result is a sparse matrix <em>B</em>.
/*! If <em>indexFrom1</em>==<b>true</b>,  <em>B</em>(<em>i,perm</em>[<em>j-1</em>])            = <em>A</em>(<em>i,j</em>) for <em>i=1,...,nr</em> and <em>j=1,...,nc</em>.
 *  If <em>indexFrom1</em>==<b>false</b>, <em>B</em>(<em>i,perm</em>[<em>j-1</em>]+<em>1</em>) = <em>A</em>(<em>i,j</em>) for <em>i=1,...,nr</em> and <em>j=1,...,nc</em>.
 *
 *  <ul>
 *  <li>  <em>A</em> is the <em>nr</em>-by-<em>nc</em> input  sparse matrix in CSC format stored in <em>store</em>[],  <em>rowIdx</em>[],  <em>colPtr</em>[].
 *  <li>  <em>B</em> is the <em>nr</em>-by-<em>nc</em> output sparse matrix in CSC format stored in <em>store2</em>[], <em>rowIdx2</em>[], <em>colPtr2</em>[].
 *  </ul>
 *
 *  \param  nc        (as input) is the number of columns of <em>A</em>.
 *  \param  perm[]    (as input) is the permutation array.
 *  \param  store[]   (as input) stores the elements of <em>A</em> in the order as <em>rowIdx</em>[].
 *  \param  rowIdx[]  (as input) stores the row indices of <em>A</em>.
 *  \param  colPtr[]  (as input) stores the column pointers of <em>A</em>.
 *  \param  store2[]  (as output) stores the elements of <em>B</em> in the order as <em>rowIdx2</em>[].
 *  \param  rowIdx2[] (as output) stores the row indices of <em>B</em>.
 *  \param  colPtr2[] (as output) stores the column pointers of <em>B</em>.
 *  \param  indexFrom1   (as input) tells whether the indices <em>perm</em>[] start with <em>1</em> or <em>0</em>. \n
 *                    If <em>indexFrom1</em>==<b>true</b>,  the indices <em>perm</em>[] range over <em>1,...,nc</em>. \n
 *                    If <em>indexFrom1</em>==<b>false</b>, the indices <em>perm</em>[] range over <em>0,...,nc-1</em>.
 *  \param  patternOnly  (as input) tells whether only the sparsity pattern is considered. \n
 *                    If <em>patternOnly</em>==<b>true</b>, then only the sparsity pattern is considered, and
 *                    <em>store</em>[] and <em>store2</em>[] will not be used.
 *
 *  \remark  The number of rows of <em>A</em>, <em>nr</em>, is not required in the computation, so it is not passed.
 *  \remark  If <em>store2</em>==<tt>NULL</tt> (<em>rowIdx2</em>==<tt>NULL</tt>, or <em>colPtr2</em>==<tt>NULL</tt>) is given, then the required memory will be allocated;
 *           otherwise it is assumed that sufficient memory has been allocated, with the address pointed to by <em>store2</em> (<em>rowIdx2</em>, <em>colPtr2</em>, respectively).
 *  \remark  This routine resembles <em>B</em>(:,<em>perm</em>)=<em>A</em> in <b>OCTAVE</b>/<b>MATLAB</b>.
 *  \remark  This routine can be used to permute rows of a sparse matrix in CSR format.
 */
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
                         Real *&store2, mkIndex *&rowIdx2, mkIndex *&colPtr2, bool indexFrom1=true, bool patternOnly=false);

//! Permute rows of a sparse matrix <em>A</em> in CSC format (not in-place). The result is a sparse matrix <em>B</em>.
/*! If <em>indexFrom1</em>==<b>true</b>,  <em>B</em>(<em>perm</em>[<em>i-1</em>],<em>j</em>)   = <em>A</em>(<em>i,j</em>) for <em>i=1,...,nr</em> and <em>j=1,...,nc</em>.
 *  If <em>indexFrom1</em>==<b>false</b>, <em>B</em>(<em>perm</em>[<em>i-1</em>]+<em>1,j</em>) = <em>A</em>(<em>i,j</em>) for <em>i=1,...,nr</em> and <em>j=1,...,nc</em>.
 *
 *  <ul>
 *  <li>  <em>A</em> is the <em>nr</em>-by-<em>nc</em> input  sparse matrix in CSC format stored in <em>store</em>[],  <em>rowIdx</em>[],  <em>colPtr</em>[].
 *  <li>  <em>B</em> is the <em>nr</em>-by-<em>nc</em> output sparse matrix in CSC format stored in <em>store2</em>[], <em>rowIdx2</em>[], <em>colPtr2</em>[].
 *  </ul>
 *
 *  \param  nc        (as input) is the number of columns of <em>A</em>.
 *  \param  perm[]    (as input) is the permutation array.
 *  \param  store[]   (as input) stores the elements of <em>A</em> in the order as <em>rowIdx</em>[].
 *  \param  rowIdx[]  (as input) stores the row indices of <em>A</em>.
 *  \param  colPtr[]  (as input) stores the column pointers of <em>A</em>.
 *  \param  store2[]  (as output) stores the elements of <em>B</em> in the order as <em>rowIdx2</em>[].
 *  \param  rowIdx2[] (as output) stores the row indices of <em>B</em>.
 *  \param  colPtr2[] (as output) stores the column pointers of <em>B</em>.
 *  \param  indexFrom1   (as input) tells whether the indices <em>perm</em>[] start with <em>1</em> or <em>0</em>. \n
 *                    If <em>indexFrom1</em>==<b>true</b>,  the indices <em>perm</em>[] range over <em>1,...,nr</em>. \n
 *                    If <em>indexFrom1</em>==<b>false</b>, the indices <em>perm</em>[] range over <em>0,...,nr-1</em>.
 *  \param  patternOnly  (as input) tells whether only the sparsity pattern is considered. \n
 *                    If <em>patternOnly</em>==<b>true</b>, then only the sparsity pattern is considered, and <em>store</em>[] and <em>store2</em>[] will not be used.
 *
 *  \remark  The number of rows of <em>A</em>, <em>nr</em>, is not required in the computation, so it is not passed.
 *  \remark  If <em>store2</em>==<tt>NULL</tt> (<em>rowIdx2</em>==<tt>NULL</tt>, or <em>colPtr2</em>==<tt>NULL</tt>) is given, then the required memory will be allocated;
 *           otherwise it is assumed that sufficient memory has been allocated, with the address pointed to by <em>store2</em> (<em>rowIdx2</em>, <em>colPtr2</em>, respectively).
 *  \remark  The elements in each column of <em>B</em> may not be sorted with respect to the row indices, even if <em>A</em> has this property.
 *           If it is desired to have sorted elements in each column, then use sort_csc_elements_in_each_column() after permute_csc_rows().
 *  \remark  This routine resembles <em>B</em>(<em>perm</em>,:)=<em>A</em> in <b>OCTAVE</b>/<b>MATLAB</b>.
 *  \remark  This routine can be used to permute columns of a sparse matrix in CSR format.
 */
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
void permute_csc_rows(mkIndex nc, const mkIndex *perm, const Real *store, const mkIndex *rowIdx, const mkIndex *colPtr,
                      Real *&store2, mkIndex *&rowIdx2, mkIndex *&colPtr2, bool indexFrom1=true, bool patternOnly=false);

//! Permute rows and columns of a sparse matrix <em>A</em> in CSC format (not in-place). The result is a sparse matrix <em>B</em>.
/*! If <em>indexFrom1</em>==<b>true</b>,  <em>B</em>(<em>rperm</em>[<em>i-1</em>]           ,<em>cperm</em>[<em>j-1</em>])            = <em>A</em>(<em>i,j</em>) for <em>i=1,...,nr</em> and <em>j=1,...,nc</em>.
 *  If <em>indexFrom1</em>==<b>false</b>, <em>B</em>(<em>rperm</em>[<em>i-1</em>]+<em>1</em>,<em>cperm</em>[<em>j-1</em>]+<em>1</em>) = <em>A</em>(<em>i,j</em>) for <em>i=1,...,nr</em> and <em>j=1,...,nc</em>.
 *
 *  <ul>
 *  <li>  <em>A</em> is the <em>nr</em>-by-<em>nc</em> input  sparse matrix in CSC format stored in <em>store</em>[],  <em>rowIdx</em>[],  <em>colPtr</em>[].
 *  <li>  <em>B</em> is the <em>nr</em>-by-<em>nc</em> output sparse matrix in CSC format stored in <em>store2</em>[], <em>rowIdx2</em>[], <em>colPtr2</em>[].
 *  </ul>
 *
 *  \param  nc        (as input) is the number of columns of <em>A</em>.
 *  \param  rperm[]   (as input) is the row permutation array.
 *  \param  cperm[]   (as input) is the column permutation array.
 *  \param  store[]   (as input) stores the elements of <em>A</em> in the order as <em>rowIdx</em>[].
 *  \param  rowIdx[]  (as input) stores the row indices of <em>A</em>.
 *  \param  colPtr[]  (as input) stores the column pointers of <em>A</em>.
 *  \param  store2[]  (as output) stores the elements of <em>B</em> in the order as <em>rowIdx2</em>[].
 *  \param  rowIdx2[] (as output) stores the row indices of <em>B</em>.
 *  \param  colPtr2[] (as output) stores the column pointers of <em>B</em>.
 *  \param  indexFrom1   (as input) tells whether the indices <em>rperm</em>[] and <em>cperm</em>[] start with <em>1</em> or <em>0</em>. \n
 *                    If <em>indexFrom1</em>==<b>true</b>,  the indices <em>rperm</em>[] and <em>cperm</em>[] range over <em>1,...,nc</em>. \n
 *                    If <em>indexFrom1</em>==<b>false</b>, the indices <em>rperm</em>[] and <em>cperm</em>[] range over <em>0,...,nc-1</em>.
 *  \param  patternOnly  (as input) tells whether only the sparsity pattern is considered. \n
 *                    If <em>patternOnly</em>==<b>true</b>, then only the sparsity pattern is considered, and
 *                    <em>store</em>[] and <em>store2</em>[] will not be used.
 *
 *  \remark  The number of rows of <em>A</em>, <em>nr</em>, is not required in the computation, so it is not passed.
 *  \remark  If <em>store2</em>==<tt>NULL</tt> (<em>rowIdx2</em>==<tt>NULL</tt>, or <em>colPtr2</em>==<tt>NULL</tt>) is given, then the required memory will be allocated;
 *           otherwise it is assumed that sufficient memory has been allocated, with the address pointed to by <em>store2</em> (<em>rowIdx2</em>, <em>colPtr2</em>, respectively).
 *  \remark  The elements in each column of <em>B</em> may not be sorted with respect to the row indices, even if <em>A</em> has this property.
 *           If it is desired to have sorted elements in each column, then use sort_csc_elements_in_each_column() after permute_csc_rows().
 *  \remark  This routine resembles <em>B</em>(<em>rperm</em>,<em>cperm</em>)=<em>A</em> in <b>OCTAVE</b>/<b>MATLAB</b>.
 *  \remark  This routine can be used to permute columns and rows of a sparse matrix in CSR format.
 */
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
                                  Real *&store2, mkIndex *&rowIdx2, mkIndex *&colPtr2, bool indexFrom1=true, bool patternOnly=false);

//! Form a tridiagonal matrix <em>A</em> as a sparse matrix with three diagonals <em>ldiag</em>[], <em>diag</em>[], and <em>udiag</em>[].
/*! The lower subdiagonal, the main diagonal, and the upper subdiagonal of <em>A</em> are formed by
 *  <em>ldiag</em>[<em>0,...,nrc-2</em>], <em>diag</em>[<em>0,...,nrc-1</em>], and <em>udiag</em>[<em>0,...,nrc-2</em>], respectively.
 *
 *  To be precise,
 *  <ul>
 *  <li>  <em>A</em>(<em>i+1,i</em>) = <em>ldiag</em>[<em>i-1</em>] for <em>i=1,...,nrc-1</em>.
 *  <li>  <em>A</em>(<em>i,i</em>)   =  <em>diag</em>[<em>i-1</em>] for <em>i=1,...,nrc</em>.
 *  <li>  <em>A</em>(<em>i,i+1</em>) = <em>udiag</em>[<em>i</em>]   for <em>i=1,...,nrc-1</em>.
 *  </ul>
 *
 *  The output sparse matrix <em>A</em> is stored in CSC format in <em>sa</em>[], <em>rowIdx</em>[], <em>colPtr</em>[].
 *
 *  \param  nrc       (as input) is the number of rows/columns of <em>A</em>.
 *  \param  ldiag[]   (as input) gives the lower subdiagonl elements of <em>A</em>.
 *  \param  diag[]    (as input) gives the main diagonal elements of <em>A</em>.
 *  \param  udiag[]   (as input) gives the upper subdiagonal elements of <em>A</em>.
 *  \param  sa[]      (as output) stores the elements of <em>A</em> in the order as <em>rowIdx</em>[].
 *  \param  rowIdx[]  (as output) stores the row indices of <em>A</em>.
 *  \param  colPtr[]  (as output) stores the column pointers of <em>A</em>.
 *
 *  \remark  If <em>sa</em>==<tt>NULL</tt> (<em>rowIdx</em>==<tt>NULL</tt>, or <em>colPtr</em>==<tt>NULL</tt>) is given, then the required memory will be allocated;
 *           otherwise it is assumed that sufficient memory has been allocated, with the address pointed to by <em>sa</em> (<em>rowIdx</em>, <em>colPtr</em>, respectively).
 */
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
                             Real *&sa, mkIndex *&rowIdx, mkIndex *&colPtr);


#ifdef USE_MEX
// convert SparseMatrix to mxArray and vice versa
void mxArray_to_spmatrix(const mxArray *mA, mkIndex nr, mkIndex nc, Real *&sa, mkIndex *&rowIdx, mkIndex *&colPtr);
mxArray *spmatrix_to_mxArray(const Real *sa, const mkIndex *rowIdx, const mkIndex *colPtr);
#endif

/** @} */  // end of group_spmatrix



////////////////////////////////////////////////////////////////////////////
//    vector - vector operations
////////////////////////////////////////////////////////////////////////////

/** @defgroup group_vector_vector Vector - Vector Operations
 *
 *  This module contains routines `without' classes related to vector - vector operations
 *  (declared in matkitfunc.h).
 *
 *  @{
 */

//! Concatenate two vectors <em>v</em> and <em>w</em>. The result is a vector <em>u</em>.
/*! To be precise, this routine computes <em>u</em>(<em>i</em>)=<em>v</em>(<em>i</em>) for <em>i=1,...,nv</em> and
 *                                  <em>u</em>(<em>nv+j</em>)<em>=w</em>(<em>j</em>) for <em>j=1,...,nw</em>.
 *  <ul>
 *  <li>  <em>v, w</em> are the input vectors of lengths <em>nv, nw</em> stored in <em>sv</em>[], <em>sw</em>[], respectively.
 *  <li>  <em>u</em>   is the output vector of length <em>nv+nw</em> stored in <em>su</em>[].
 *  </ul>
 *
 *  \param  nv        (as input) is the length of <em>v</em>.
 *  \param  sv[]      (as input) stores the elements of <em>v</em>.
 *  \param  nw        (as input) is the length of <em>w</em>.
 *  \param  sw[]      (as input) stores the elements of <em>w</em>.
 *  \param  su[]      (as output) stores the elements of <em>u</em>.
 *
 *  \remark  If <em>su</em>==<tt>NULL</tt> is given, then the required memory of size <em>nv+nw</em> will be allocated for the output vector <em>u</em>;
 *           otherwise it is assumed that sufficient memory has been allocated, with the address pointed to by <em>su</em>.
 */
// concatenate two vectors v and w to form u, where
//    v, w are the input vectors of lengths nv,nw stored in sv[], sw[], respectively, and
//    u is the output vector u length nv+nw stored in su[]
// to be precise, this routine computes u(i)=v(i) for i=1,...,nv and u(nv+j)=w(j) for j=1,...,nw
// if su==NULL is given, then the required memory of size nv+nw will be allocated for the output vector u;
//    otherwise it is assumed that sufficient memory has been allocated, with the address pointed to by su
void vector_concatenate(mkIndex nv, const Real *sv, mkIndex nw, const Real *sw, Real *&su);

//! Return the inner product <em>v'*w</em> of two vectors <em>v</em> and <em>w</em>.
/*! <ul>
 *  <li>  <em>v, w</em> are the input vectors of length <em>len</em> stored in <em>sv</em>[], <em>sw</em>[], respectively.
 *  </ul>
 *
 *  \param  len       (as input) is the length of <em>v, w</em>.
 *  \param  sv[]      (as input) stores the elements of <em>v</em>.
 *  \param  sw[]      (as input) stores the elements of <em>w</em>.
 *
 *  \return  The value of <em>v</em>(<em>1</em>)*<em>w</em>(<em>1</em>)<em>+...+v</em>(<em>len</em>)<em>*w</em>(<em>len</em>).
 *  \remark  This routine invokes the level 1 <b>BLAS</b> routine <b>xDOT</b> if <tt>USE_BLAS</tt> is defined.
 */
// return the inner product v'*w of two vectors v and w, where
//    v, w are the input vectors of length len stored in sv[], sw[], respectively
// to be precise, the return value is v(1)*w(1)+...+v(len)*w(len)
// this routine invokes the level 1 BLAS routine xDOT if USE_BLAS is defined
Real vector_inner_product(mkIndex len, const Real *sv, const Real *sw);

//! Compute <em>A=alp*v*w'</em> (if <em>init0</em>==<b>true</b> or <em>sa</em>==<tt>NULL</tt>) or <em>A+=alp*v*w'</em> (if <em>init0</em>==<b>false</b> and <em>sa</em>\f$\neq\f$<tt>NULL</tt>), where <em>A</em> is a general matrix, <em>v, w</em> are vectors, and <em>alp</em> is a scalar.
/*! <ul>
 *  <li>  <em>v, w</em> are the input vectors of lengths <em>nr, nc</em> stored in <em>sv</em>[], <em>sw</em>[], respectively.
 *  <li>  <em>A</em> is the <em>nr</em>-by-<em>nc</em> output general matrix stored in <em>sa</em>[].
 *  </ul>
 *
 *  \param  nr        (as input) is length of <em>v</em>, which is also the number of rows of <em>A</em>.
 *  \param  nc        (as input) is length of <em>w</em>, which is also the number of columns of <em>A</em>.
 *  \param  alp       (as input) is a scalar.
 *  \param  sv[]      (as input) stores the elements of v.
 *  \param  sw[]      (as input) stores the elements of w.
 *  \param  sa[]      (as output) stores the elements of A.
 *  \param  init0     (as input) tells whether to initialize <em>A</em> with zeros or not.
 *                    If <em>sa</em>==<tt>NULL</tt> is given, then <em>A</em> will be initialized with zeros anyway.
 *
 *  \remark  If <em>sa</em>==<tt>NULL</tt> is given, then the required memory of size <em>nr*nc</em> will be allocated for the output matrix <em>A</em>;
 *           otherwise it is assumed that sufficient memory has been allocated, with the address pointed to by <em>sa</em>.
 *  \remark  This routine invokes the level 2 <b>BLAS</b> routine <b>xGER</b> if <tt>USE_BLAS</tt> is defined.
 */
// compute A = alp*v*w' (if init0==true or sa==NULL) or A += alp*v*w' (if init0==false and sa!=NULL), where
//    alp is a scalar,
//    v, w are the input vectors of lengths nr, nc stored in sv[], sw[], respectively, and
//    A is the nr-by-nc output general matrix stored in sa[]
// if sa==NULL is given, then the required memory of size nr*nc will be allocated for the output matrix A;
//    otherwise it is assumed that sufficient memory has been allocated, with the address pointed to by sa
// this routine invokes the level 2 BLAS routine xGER if USE_BLAS is defined
void vector_outer_product(mkIndex nr, mkIndex nc, Real alp, const Real *sv, const Real *sw, Real *&sa, bool init0=true);

//! Compute <em>A=alp*v*v'</em> (if <em>init0</em>==<b>true</b> or <em>sa</em>==<tt>NULL</tt>) or <em>A+=alp*v*v'</em> (if <em>init0</em>==<b>false</b> and <em>sa</em>\f$\neq\f$<tt>NULL</tt>), where <em>A</em> is a symmetric matrix, <em>v</em> is a vector, and <em>alp</em> is a scalar.
/*! <ul>
 *  <li>  <em>v</em> is the input vector of length <em>nrc</em> stored in <em>sv</em>[].
 *  <li>  <em>A</em> is the <em>nrc</em>-by-<em>nrc</em> output symmetric matrix stored in packed form in <em>sa</em>[].
 *  </ul>
 *
 *  \param  nrc       (as input) is the length of <em>v</em>, which is also the number of rows/columns of <em>A</em>.
 *  \param  alp       (as input) is a scalar.
 *  \param  sv[]      (as input) stores the elements of <em>v</em>.
 *  \param  sa[]      (as output) stores the elements of <em>A</em>.
 *  \param  init0     (as input) tells whether to initialize <em>A</em> with zeros or not.
 *                    If <em>sa</em>==<tt>NULL</tt> is given, then <em>A</em> will be initialized with zeros anyway.
 *
 *  \remark  If <em>sa</em>==<tt>NULL</tt> is given, then the required memory of size <em>nrc*(nrc+1)/2</em> will be allocated for the output symmetric matrix <em>A</em>;
 *           otherwise it is assumed that sufficient memory has been allocated, with the address pointed to by <em>sa</em>.
 *  \remark  This routine invokes the level 2 <b>BLAS</b> routine <b>xSPR</b> if <tt>USE_BLAS</tt> is defined.
 */
// compute A = alp*v*v' (if init0==true or sa==NULL) or A += alp*v*v' (if init0==false and sa!=NULL), where
//    alp is a scalar,
//    v is the input vector of length nrc stored in sv[], and
//    A is the nrc-by-nrc output symmetric matrix stored in packed form in sa[]
// if sa==NULL is given, then the required memory of size nrc*(nrc+1)/2 for will be allocated the output symmetric matrix A;
//    otherwise it is assumed that sufficient memory has been allocated, with the address pointed to by sa
// this routine invokes the level 2 BLAS routine xSPR if USE_BLAS is defined
void vector_outer_product(mkIndex nrc, Real alp, const Real *sv, Real *&sa, bool init0=true);

//! Compute <em>A=alp*(v*w'+w*v')</em> (if <em>init0</em>==<b>true</b> or <em>sa</em>==<tt>NULL</tt>) or <em>A+=alp*(v*w'+w*v')</em> (if <em>init0</em>==<b>false</b> and <em>sa</em>\f$\neq\f$<tt>NULL</tt>), where <em>A</em> is a symmetric matrix, <em>v, w</em> are vectors, and <em>alp</em> is a scalar.
/*! <ul>
 *  <li>  <em>v, w</em> are the input vectors of length <em>nrc</em> stored in <em>sv</em>[], <em>sw</em>[], respectively.
 *  <li>  <em>A</em> is the <em>nrc</em>-by-<em>nrc</em> output symmetric matrix stored in packed form in <em>sa</em>[].
 *  </ul>
 *
 *  \param  nrc       (as input) is length of <em>v, w</em>, which is also the number of rows/columns of <em>A</em>.
 *  \param  alp       (as input) is a scalar.
 *  \param  sv[]      (as input) stores the elements of v.
 *  \param  sw[]      (as input) stores the elements of w.
 *  \param  sa[]      (as output) stores the elements of A.
 *  \param  init0     (as input) tells whether to initialize <em>A</em> with zeros or not.
 *                    If <em>sa</em>==<tt>NULL</tt> is given, then <em>A</em> will be initialized with zeros anyway.
 *
 *  \remark  If <em>sa</em>==<tt>NULL</tt> is given, then the required memory of size <em>nrc*(nrc+1)/2</em> will be allocated for the output matrix <em>A</em>;
 *           otherwise it is assumed that sufficient memory has been allocated, with the address pointed to by <em>sa</em>.
 *  \remark  This routine invokes the level 2 <b>BLAS</b> routine <b>xSPR2</b> if <tt>USE_BLAS</tt> is defined.
 */
// compute A = alp*(v*w'+w*v') (if init0==true or sa==NULL) or A += alp*(v*w'+w*v') (if init0==false and sa!=NULL), where
//    alp is a scalar, and
//    v, w are the input vectors of length nrc stored in sv[], sw[], respectively
//    A is the nrc-by-nrc output symmetric matrix stored in packed form in sa[]
// if sa==NULL is given, then the required memory of size len*(len+1)/2 will be allocated for the output symmetric matrix A;
//    otherwise it is assumed that sufficient memory has been allocated, with the address pointed to by sa
// this routine invokes the level 2 BLAS routine xSPR2 if USE_BLAS is defined
void vector_symmetric_outer_product(mkIndex nrc, Real alp, const Real *sv, const Real *sw, Real *&sa, bool init0=true);

/** @} */  // end of group_vector_vector



////////////////////////////////////////////////////////////////////////////
//    matrix - vector operations
////////////////////////////////////////////////////////////////////////////

/** @defgroup group_matrix_vector Matrix - Vector Operations
 *
 *  This module contains routines `without' classes related to matrix - vector operations
 *  (declared in matkitfunc.h).
 *
 *  @{
 */

//! Compute <em>w=A*v</em> (if <em>init0</em>==<b>true</b> or <em>sw</em>==<tt>NULL</tt>) or <em>w+=A*v</em> (if <em>init0</em>==<b>false</b> and <em>sw</em>\f$\neq\f$<tt>NULL</tt>), where <em>A</em> is a general matrix and <em>v, w</em> are vectors.
/*! <ul>
 *  <li>  <em>A</em> is the <em>nr</em>-by-<em>nc</em> input general matrix stored in <em>sa</em>[].
 *  <li>  <em>v</em> is the input vector of length <em>nc</em> stored in <em>sv</em>[].
 *  <li>  <em>w</em> is the output vector of length <em>nr</em> stored in <em>sw</em>[].
 *  </ul>
 *
 *  \param  nr        (as input) is the number of rows of <em>A</em>, which is also the length of <em>w</em>.
 *  \param  nc        (as input) is the number of columns of <em>A</em>, which is also the length of <em>v</em>.
 *  \param  sa[]      (as input) stores the elements of <em>A</em>.
 *  \param  sv[]      (as input) stores the elements of <em>v</em>.
 *  \param  sw[]      (as output) stores the elements of <em>w</em>.
 *  \param  init0     (as input) tells whether to initialize <em>w</em> with zeros or not.
 *                    If <em>sw</em>==<tt>NULL</tt> is given, then <em>w</em> will be initialized with zeros anyway.
 *
 *  \remark  If <em>sw</em>==<tt>NULL</tt> is given, then the required memory of size <em>nr</em> will be allocated for the output vector <em>w</em>;
 *           otherwise, it is assumed that sufficient memory has been allocated, with the address pointed to by <em>sw</em>.
 *  \remark  This routine invokes the level 2 <b>BLAS</b> routine <b>xGEMV</b> if <tt>USE_BLAS</tt> is defined.
 */
// compute w = A*v (if init0==true or sw==NULL) or w += A*v (if init0==false and sw!=NULL), where
//    A is the nr-by-nc input general matrix stored in sa[],
//    v is the input vector of length nc stored in sv[], and
//    w is the output vector of length nr stored in sw[]
// if sw==NULL is given, then the required memory of size nr will be allocated for the output vector w;
//    otherwise it is assumed that sufficient memory has been allocated, with the address pointed to by sw
// this routine invokes the level 2 BLAS routine xGEMV if USE_BLAS is defined
void matrix_vector_multiplication(mkIndex nr, mkIndex nc, const Real *sa, const Real *sv, Real *&sw, bool init0=true);

//! Compute <em>w'=v'*A</em> (if <em>init0</em>==<b>true</b> or <em>sw</em>==<tt>NULL</tt>) or <em>w'+=v'*A</em> (if <em>init0</em>==<b>false</b> and <em>sw</em>\f$\neq\f$<tt>NULL</tt>), where <em>A</em> is a general matrix and <em>v, w</em> are vectors.
/*! <ul>
 *  <li>  <em>A</em> is the <em>nr</em>-by-<em>nc</em> input general matrix stored in <em>sa</em>[].
 *  <li>  <em>v</em> is the input vector of length <em>nr</em> stored in <em>sv</em>[].
 *  <li>  <em>w</em> is the output vector of length <em>nc</em> stored in <em>sw</em>[].
 *  </ul>
 *
 *  \param  nr        (as input) is the number of rows of <em>A</em>, which is also the length of <em>v</em>.
 *  \param  nc        (as input) is the number of columns of <em>A</em>, which is also the length of <em>w</em>.
 *  \param  sa[]      (as input) stores the elements of <em>A</em>.
 *  \param  sv[]      (as input) stores the elements of <em>v</em>.
 *  \param  sw[]      (as output) stores the elements of <em>w</em>.
 *  \param  init0     (as input) tells whether to initialize <em>w</em> with zeros or not.
 *                    If <em>sw</em>==<tt>NULL</tt> is given, then <em>w</em> will be initialized with zeros anyway.
 *
 *  \remark  If <em>sw</em>==<tt>NULL</tt> is given, then the required memory of size <em>nc</em> will be allocated for the output vector <em>w</em>;
 *           otherwise, it is assumed that sufficient memory has been allocated, with the address pointed to by <em>sw</em>.
 *  \remark  This routine invokes the level 2 <b>BLAS</b> routine <b>xGEMV</b> if <tt>USE_BLAS</tt> is defined.
 */
// compute w' = v'*A (if init0==true or sw==NULL) or w' += v'*A (if init0==false and sw!=NULL), where
//    A is the nr-by-nc input general matrix stored in sa[],
//    v is the input vector of length nr stored in sv[], and
//    w is the output vector of length nc stored in sw[]
// if sw==NULL is given, then the required memory of size nc will be allocated for the output vector w;
//    otherwise it is assumed that sufficient memory has been allocated, with the address pointed to by sw
// this routine invokes the level 2 BLAS routine xGEMV if USE_BLAS is defined
void vector_matrix_multiplication(mkIndex nr, mkIndex nc, const Real *sv, const Real *sa, Real *&sw, bool init0=true);

/** @} */  // end of group_matrix_vector


////////////////////////////////////////////////////////////////////////////
//    symmetric matrix - vector operations
////////////////////////////////////////////////////////////////////////////

/** @defgroup group_symmatrix_vector Symmetric Matrix - Vector Operations
 *
 *  This module contains routines `without' classes related to symmetric matrix - vector operations
 *  (declared in matkitfunc.h).
 *
 *  @{
 */

//! Compute <em>w=A*v</em> (if <em>init0</em>==<b>true</b> or <em>sw</em>==<tt>NULL</tt>) or <em>w+=A*v</em> (if <em>init0</em>==<b>false</b> and <em>sw</em>\f$\neq\f$<tt>NULL</tt>), where <em>A</em> is a symmetric matrix and <em>v, w</em> are vectors.
/*! <ul>
 *  <li>  <em>A</em> is the <em>nrc</em>-by-<em>nrc</em> symmetric matrix stored in packed form in <em>sa</em>[].
 *  <li>  <em>v</em> is the input vector of length <em>nrc</em> stored in <em>sv</em>[].
 *  <li>  <em>w</em> is the output vector of length <em>nrc</em> stored in <em>sw</em>[].
 *  </ul>
 *
 *  \param  nrc       (as input) is the number of rows/columns of <em>A</em>, which is also the length of <em>v, w</em>.
 *  \param  sa[]      (as input) stores the elements of <em>A</em>.
 *  \param  sv[]      (as input) stores the elements of <em>v</em>.
 *  \param  sw[]      (as output) stores the elements of <em>w</em>.
 *  \param  init0     (as input) tells whether to initialize <em>w</em> with zeros or not.
 *                    If <em>sw</em>==<tt>NULL</tt> is given, then <em>w</em> will be initialized with zeros anyway.
 *
 *  \remark  If <em>sw</em>==<tt>NULL</tt> is given, then the required memory of size <em>nrc</em> will be allocated for the output vector <em>w</em>;
 *           otherwise, it is assumed that sufficient memory has been allocated, with the address pointed to by <em>sw</em>.
 *  \remark  This routine invokes the level 2 <b>BLAS</b> routine <b>xSPMV</b> if <tt>USE_BLAS</tt> is defined
 */
// compute w = A*v (if init0==true or sw==NULL) or w' += v'*A (if init0==false and sw!=NULL), where
//    A is the nrc-by-nrc input symmetric matrix stored in packed form in sa[],
//    v is the input vector of length nrc stored in sv[], and
//    w is the output vector of length nrc stored in sw[]
// if sw==NULL is given, then the required memory of size nrc will be allocated for the output vector w;
//    otherwise it is assumed that sufficient memory has been allocated, with the address pointed to by sw
// this routine invokes the level 2 BLAS routine xSPMV if USE_BLAS is defined
void symmatrix_vector_multiplication(mkIndex nrc, const Real *sa, const Real *sv, Real *&sw, bool init0=true);

/** @} */  // end of group_symmatrix_vector



////////////////////////////////////////////////////////////////////////////
//    sparse matrix - vector operations
////////////////////////////////////////////////////////////////////////////

/** @defgroup group_spmatrix_vector Sparse Matrix - Vector Operations
 *
 *  This module contains routines `without' classes related to sparse matrix - vector operations
 *  (declared in matkitfunc.h).
 *
 *  @{
 */

//! Compute <em>w=A*v</em> (if <em>init0</em>==<b>true</b> or <em>sw</em>==<tt>NULL</tt>) or <em>w+=A*v</em> (if <em>init0</em>==<b>false</b> and <em>sw</em>\f$\neq\f$<tt>NULL</tt>), where <em>A</em> is a sparse matrix and <em>v, w</em> are vectors.
/*! <ul>
 *  <li>  <em>A</em> is the <em>nr</em>-by-<em>nc</em> input sparse matrix in CSC format stored in <em>store</em>[], <em>rowIdx</em>[], <em>colPtr</em>[].
 *  <li>  <em>v</em> is the input vector of length <em>nc</em> stored in <em>sv</em>[].
 *  <li>  <em>w</em> is the output vector of length <em>nr</em> stored in <em>sw</em>[].
 *  </ul>
 *
 *  \param  nr        (as input) is the number of rows of <em>A</em>, which is also the length of <em>w</em>.
 *  \param  nc        (as input) is the number of columns of <em>A</em>, which is also the length of <em>v</em>.
 *  \param  store[]   (as input) stores the elements of <em>A</em>.
 *  \param  rowIdx[]  (as input) stores the row indices of <em>A</em>.
 *  \param  colPtr[]  (as input) stores the column pointers of <em>A</em>.
 *  \param  sv[]      (as input) stores the elements of <em>v</em>.
 *  \param  sw[]      (as output) stores the elements of <em>w</em>.
 *  \param  init0     (as input) tells whether to initialize <em>w</em> with zeros or not.
 *                    If <em>sw</em>==<tt>NULL</tt> is given, then <em>w</em> will be initialized with zeros anyway.
 *
 *  \remark  If <em>sw</em>==<tt>NULL</tt> is given, then the required memory of size <em>nr</em> will be allocated for the output vector <em>w</em>;
 *           otherwise, it is assumed that sufficient memory has been allocated, with the address pointed to by <em>sw</em>.
 *  \remark  This routine can be used to compute <em>w'=v'*A</em> or <em>w'+=v'*A</em> with <em>A</em> stored in CSR format.
 */
// compute w = A*v (if init0==true or sw==NULL) or w += A*v (if init0==false and sw!=NULL), where
//    A is the nr-by-nc input sparse matrix in CSC format stored in store[], rowIdx[], colPtr[],
//    v is the input vector of length nc stored in sv[], and
//    w is the output vector of length nr stored in sw[]
// if sw==NULL is given, then the required memory of size nr will be allocated for the output vector w;
//    otherwise it is assumed that sufficient memory has been allocated, with the address pointed to by sw
// this routine can be used to compute w' = v'*A  or w' += v'*A, where A is a sparse matrix stored in CSR format
void spmatrix_vector_multiplication(mkIndex nr, mkIndex nc, const Real *store, const mkIndex *rowIdx, const mkIndex *colPtr, const Real *sv, Real *&sw, bool init0=true);

//! Compute <em>w'=v'*A</em> (if <em>init0</em>==<b>true</b> or <em>sw</em>==<tt>NULL</tt>) or <em>w'+=v'*A</em> (if <em>init0</em>==<b>false</b> and <em>sw</em>\f$\neq\f$<tt>NULL</tt>), where <em>A</em> is a sparse matrix and <em>v, w</em> are vectors.
/*! <ul>
 *  <li>  <em>A</em> is the <em>nr</em>-by-<em>nc</em> input sparse matrix in CSC format stored in <em>store</em>[], <em>rowIdx</em>[], <em>colPtr</em>[].
 *  <li>  <em>v</em> is the input vector of length <em>nr</em> stored in <em>sv</em>[].
 *  <li>  <em>w</em> is the output vector of length <em>nc</em> stored in <em>sw</em>[].
 *  </ul>
 *
 *  \param  nc        (as input) is the number of columns of <em>A</em>, which is also the length of <em>w</em>.
 *  \param  store[]   (as input) stores the elements of <em>A</em>.
 *  \param  rowIdx[]  (as input) stores the row indices of <em>A</em>.
 *  \param  colPtr[]  (as input) stores the column pointers of <em>A</em>.
 *  \param  sv[]      (as input) stores the elements of <em>v</em>.
 *  \param  sw[]      (as output) stores the elements of <em>w</em>.
 *  \param  init0     (as input) tells whether to initialize <em>w</em> with zeros or not.
 *                    If <em>sw</em>==<tt>NULL</tt> is given, then <em>w</em> will be initialized with zeros anyway.
 *
 *  \remark  The number of rows in <em>A</em>, <em>nr</em>, is not required in the computation, so it is not passed.
 *  \remark  If <em>sw</em>==<tt>NULL</tt> is given, then the required memory of size <em>nr</em> will be allocated for the output vector <em>w</em>;
 *           otherwise, it is assumed that sufficient memory has been allocated, with the address pointed to by <em>sw</em>.
 *  \remark  This routine can be used to compute <em>w=A*v</em> or <em>w+=A*v</em>, where <em>A</em> is a sparse matrix stored in CSR format.
 */
// compute w' = v'*A (if init0==true or sw==NULL) or w' += v'*A (if init0==false and sw!=NULL), where
//    A is the nr-by-nc input sparse matrix in CSC format stored in store[], rowIdx[], colPtr[],
//    v is the input vector of length nr stored in sv[], and
//    w is the output vector of length nc stored in sw[]
// note that the number of rows in A, nr, is not required in the computation, so it is not passed
// if sw==NULL is given, then the required memory of size nc will be allocated for the output vector w;
//    otherwise it is assumed that sufficient memory has been allocated, with the address pointed to by sw
// this routine can be used to compute w = A*v  or w += A*v, where A is a sparse matrix stored in CSR format
void vector_spmatrix_multiplication(mkIndex nc, const Real *sv, const Real *store, const mkIndex *rowIdx, const mkIndex *colPtr, Real *&sw, bool init0=true);

/** @} */  // end of group_spmatrix_vector



////////////////////////////////////////////////////////////////////////////
//    matrix - matrix operations
////////////////////////////////////////////////////////////////////////////

/** @defgroup group_matrix_matrix  Matrix - Matrix Operations
 *
 *  This module contains routines `without' classes related to matrix - matrix operations
 *  (declared in matkitfunc.h).
 *
 *  @{
 */

//! Horizontal concatenation of two general matrices <em>A</em> and <em>B</em>. The result is a general matrix <em>C</em>.
/*! <ul>
 *  <li>  <em>A, B</em> are the <em>nr</em>-by-<em>nca</em>, <em>nr</em>-by-<em>ncb</em> input general matrices stored in <em>sa</em>[], <em>sb</em>[], respectively.
 *  <li>  <em>C</em> is the <em>nr</em>-by-(<em>nca+ncb</em>) output general matrix stored in <em>sc</em>[].
 *  </ul>
 *
 *  \param  nr        (as input) is the number of rows of <em>A, B, C</em>.
 *  \param  nca       (as input) is the number of columns of <em>A</em>.
 *  \param  sa[]      (as input) stores the elements of <em>A</em>.
 *  \param  ncb       (as input) is the number of columns of <em>B</em>.
 *  \param  sb[]      (as input) stores the elements of <em>B</em>.
 *  \param  sc[]      (as output) stores the elements of <em>C</em>.
 *
 *  \remark  If <em>sc</em>==<tt>NULL</tt> is given, then the required memory of size <em>nc*(nca+ncb)</em> will be allocated for the output matrix <em>C</em>;
 *           otherwise, it is assumed that sufficient memory has been allocated, with the address pointed to by <em>sc</em>.
 *  \remark  This static function resembles [<em>A,B</em>] or equivalently <b>horzcat</b>(<em>A,B</em>) in <b>OCTAVE</b>/<b>MATLAB</b>.
 */
// horizontal concatenation of two matrices A and B, where
//    A, B are the nr-by-nca, nr-by-ncb input general matrices stored in sa[], sb[], respectively
// A and B must have the same number of rows
// the result is an nr-by-(nca+ncb) general matrix C stored in sc[]
// if sc==NULL is given, then the required memory of size nr*(nca+ncb) will be allocated for the output matrix C;
//    otherwise it is assumed that sufficient memory has been allocated, with the address pointed to by sc
// this static function resembles [A, B] or equivalently horzcat(A,B) in OCTAVE/MATLAB
void horizontal_concatenate(mkIndex nr, mkIndex nca, const Real *sa, mkIndex ncb, const Real *sb, Real *&sc);

//! Vertical concatenation of two matrices <em>A</em> and <em>B</em>. The result is a general matrix <em>C</em>.
/*! <ul>
 *  <li>  <em>A, B</em> are the <em>nra</em>-by-<em>nc</em>, <em>nrb</em>-by-<em>nc</em> input general matrices stored in <em>sa</em>[], <em>sb</em>[], respectively.
 *  <li>  <em>C</em> is the (<em>nra+nrb</em>)-by-<em>nc</em> output general matrix stored in <em>sc</em>[].
 *  </ul>
 *
 *  \param  nra       (as input) is the number of rows of <em>A</em>.
 *  \param  nc        (as input) is the number of columns of <em>A, B, C</em>.
 *  \param  sa[]      (as input) stores the elements of <em>A</em>.
 *  \param  nrb       (as input) is the number of rows of <em>B</em>.
 *  \param  sb[]      (as input) stores the elements of <em>B</em>.
 *  \param  sc[]      (as output) stores the elements of <em>C</em>.
 *
 *  \remark  If <em>sc</em>==<tt>NULL</tt> is given, then the required memory of size (<em>nra+nrb</em>)*<em>nc</em> will be allocated for the output matrix <em>C</em>;
 *           otherwise, it is assumed that sufficient memory has been allocated, with the address pointed to by <em>sc</em>.
 *  \remark  This static function resembles [<em>A;B</em>] or equivalently <b>vertcat</b>(<em>A,B</em>) in <b>OCTAVE</b>/<b>MATLAB</b>.
 */
// vertical concatenation of two matrices A and B, where
//    A, B are the nra-by-nc, nrb-by-nc input general matrices stored in sa[], sb[], respectively
// A and B must have the same number of columns
// the result is an (nra+nrb)-by-nc general matrix C stored in sc[]
// if sc==NULL is given, then the required memory of size (nra+nrb)*nc will be allocated for the output matrix C;
//    otherwise it is assumed that sufficient memory has been allocated, with the address pointed to by sc
// this static function resembles [A; B] or equivalently vertcat(A,B) in OCTAVE/MATLAB
void vertical_concatenate(mkIndex nra, mkIndex nc, const Real *sa, mkIndex nrb, const Real *sb, Real *&sc);

//! Sum up two general matrices <em>A,B</em>, with zeros padded if <em>A</em> and <em>B</em> are of different dimensions. The result is a general matrix <em>C</em>.
/*! <ul>
 *  <li>  <em>A, B</em> are the <em>nra</em>-by-<em>nca</em>, <em>nrb</em>-by-<em>ncb</em> input general matrices stored in <em>sa</em>[], <em>sb</em>[], respectively.
 *  <li>  <em>C</em> is the <em>nrc</em>-by-<em>ncc</em> output general matrix stored in <em>sc</em>[], where <em>nrc</em>=<b>max</b>(<em>nra,nrb</em>) and <em>ncc</em>=<b>max</b>(<em>nrb,ncb</em>).
 *  </ul>
 *
 *  \param  nra       (as input) is the number of rows of <em>A</em>.
 *  \param  nca       (as input) is the number of columns of <em>A</em>.
 *  \param  sa[]      (as input) stores the elements of <em>A</em>.
 *  \param  nrb       (as input) is the number of rows of <em>B</em>.
 *  \param  ncb       (as input) is the number of columns of <em>B</em>.
 *  \param  sb[]      (as input) stores the elements of <em>B</em>.
 *  \param  sc[]      (as output) stores the elements of <em>C</em>.
 *
 *  \remark  If <em>sc</em>==<tt>NULL</tt> is given, then the required memory of size <em>nrc*ncc</em> will be allocated for the output matrix <em>C</em>;
 *           otherwise, it is assumed that sufficient memory has been allocated, with the address pointed to by <em>sc</em>.
 */
// sum up two general matrices A,B, with zeros padded if A and B are of different dimensions, where
//    A, B are the nra-by-nca, nrb-by-ncb input general matrices stored in sa[], sb[], respectively
// the result is an nrc-by-ncc general matrix C stored in sc[], where nrc=max(nra,nrb) and ncc=max(nrb,ncb)
// if sc==NULL is given, then the required memory of size nrc*ncc will be allocated for output matrix C;
//    otherwise, it is assumed that sufficient memory has been allocated, with the address pointed to by sc
void matrix_xsum(mkIndex nra, mkIndex nca, const Real *sa, mkIndex nrb, mkIndex ncb, const Real *sb, Real *&sc);

//! Compute <em>C=A*B</em> (if <em>init0</em>==<b>true</b> or <em>sc</em>==<tt>NULL</tt>) or <em>C+=A*B</em> (if <em>init0</em>==<b>false</b> and <em>sc</em>\f$\neq\f$<tt>NULL</tt>), where <em>A,B,C</em> are general matrices.
/*! <ul>
 *  <li>  <em>A, B</em> are the <em>nr</em>-by-<em>nrc</em>, <em>nrc</em>-by-<em>nc</em> input general matrices stored in <em>sa</em>[], <em>sb</em>[], respectively.
 *  <li>  <em>C</em> is the <em>nr</em>-by-<em>nc</em> output general matrix stored in <em>sc</em>[].
 *  </ul>
 *
 *  \param  nr        (as input) is the number of rows of <em>A,C</em>.
 *  \param  nrc       (as input) is the number of columns of <em>A</em>, which is also the number of rows of <em>B</em>.
 *  \param  nc        (as input) is the number of columns of <em>B,C</em>.
 *  \param  sa[]      (as input) stores the elements of <em>A</em>.
 *  \param  sb[]      (as input) stores the elements of <em>B</em>.
 *  \param  sc[]      (as output) stores the elements of <em>C</em>.
 *  \param  init0     (as input) tells whether to initialize <em>C</em> with zeros or not.
 *                    If <em>sc</em>==<tt>NULL</tt> is given, then <em>C</em> will be initialized with zeros anyway.
 *
 *  \remark  If <em>sc</em>==<tt>NULL</tt> is given, then the required memory of size <em>nr*nc</em> will be allocated for the output matrix <em>C</em>;
 *           otherwise, it is assumed that sufficient memory has been allocated, with the address pointed to by <em>sc</em>.
 *  \remark  This routine invokes the level 3 <b>BLAS</b> routine <b>xGEMM</b> if <tt>USE_BLAS</tt> is defined.
 */
// compute C = A*B (if init0==true or sc==NULL) or C += A*B (if init0==false and sc!=NULL), where
//    A, B are the nr-by-nrc, nrc-by-nc input general matrices stored in sa[], sb[], respectively, and
//    C is the nr-by-nc output general matrix stored in sc[]
// if sc==NULL is given, then the required memory nr*nc will be allocated for the output matrix C;
//    otherwise it is assumed that sufficient memory has been allocated, with the address pointed to by sc
// this routine invokes the level 3 BLAS routine xGEMM if USE_BLAS is defined
void matrix_matrix_multiplication(mkIndex nr, mkIndex nrc, mkIndex nc, const Real *sa, const Real *sb, Real *&sc, bool init0=true);

/** @} */  // end of group_matrix_matrix



////////////////////////////////////////////////////////////////////////////
//    matrix - symmetric matrix operations
////////////////////////////////////////////////////////////////////////////

/** @defgroup group_matrix_symmatrix  Matrix - Symmetric Matrix Operations
 *
 *  This module contains routines `without' classes related to matrix - symmetric matrix operations
 *  (declared in matkitfunc.h).
 *
 *  @{
 */
// matrix - symmetric matrix operations

//! Compute <em>C=A+bet*B</em> (if <em>sa\f$\neq\f$sc</em>) or <em>A+=bet*B</em> (if <em>sa==sc</em>), where <em>A,C</em> are general matrices, <em>B</em> is a symmetric matrix, and <em>bet</em> is a scalar.
/*! <ul>
 *  <li>  <em>A</em> is the <em>nrc</em>-by-<em>nrc</em> input general matrix stored in <em>sa</em>[].
 *  <li>  <em>C</em> is the <em>nrc</em>-by-<em>nrc</em> output general matrix stored in <em>sc</em>[].
 *  <li>  <em>B</em> is the <em>nrc</em>-by-<em>nrc</em> input symmetric matrix stored in packed form in <em>sb</em>[].
 *  </ul>
 *
 *  \param  nrc       (as input) is the number of rows/columns of <em>A,B,C</em>.
 *  \param  sa[]      (as input) stores the elements of <em>A</em>.
 *  \param  bet       (as input) is a scalar.
 *  \param  sb[]      (as input) stores the elements of <em>B</em>.
 *  \param  sc[]      (as output) stores the elements of <em>C</em>.
 *
 *  \remark  If <em>sc</em>==<tt>NULL</tt> is given, then the required memory of size <em>nrc*nrc</em> will be allocated for the output matrix <em>C</em>;
 *           otherwise, it is assumed that sufficient memory has been allocated, with the address pointed to by <em>sc</em>.
 */
// compute C = A + bet*B (if sa!=sc) or A += bet*B (if sa==sc), where
//    A is the nrc-by-nrc input general matrix stored in sa[],
//    B is the nrc-by-nrc input symmetric matrix stored in packed form in sb[], and
//    C is the nrc-by-nrc output general matrices stored in sc[]
// if sc==NULL is given, then the required memory of size nrc*nrc will be allocated for the output matrix C;
//    otherwise it is assumed that sufficient memory has been allocated, with the address pointed to by sc
void matrix_symmatrix_addition(mkIndex nrc, const Real *sa, Real bet, const Real *sb, Real *&sc);

//! Compute <em>C=A*B</em> (if <em>init0</em>==<b>true</b> or <em>sc</em>==<tt>NULL</tt>) or <em>C+=A*B</em> (if <em>init0</em>==<b>false</b> and <em>sc</em>\f$\neq\f$<tt>NULL</tt>), where <em>A,C</em> are general matrices and <em>B</em> is a symmetric matrix.
/*! <ul>
 *  <li>  <em>A</em> is the <em>nr</em>-by-<em>nrc</em> input general matrix stored in <em>sa</em>[].
 *  <li>  <em>B</em> is the <em>nrc</em>-by-<em>nrc</em> input symmetric matrix stored in packed form in <em>sb</em>[].
 *  <li>  <em>C</em> is the <em>nr</em>-by-<em>nrc</em> output general matrix stored in <em>sc</em>[].
 *  </ul>
 *
 *  \param  nr        (as input) is the number of rows of <em>A,C</em>.
 *  \param  nrc       (as input) is the number of columns of <em>A,C</em>, which is also the number of rows/columns of <em>B</em>.
 *  \param  sa[]      (as input) stores the elements of <em>A</em>.
 *  \param  sb[]      (as input) stores the elements of <em>B</em>.
 *  \param  sc[]      (as output) stores the elements of <em>C</em>.
 *  \param  init0     (as input) tells whether to initialize <em>C</em> with zeros or not.
 *                    If <em>sc</em>==<tt>NULL</tt> is given, then <em>C</em> will be initialized with zeros anyway.
 *
 *  \remark  If <em>sc</em>==<tt>NULL</tt> is given, then the required memory of size <em>nr*nrc</em> will be allocated for the output matrix <em>C</em>;
 *           otherwise, it is assumed that sufficient memory has been allocated, with the address pointed to by <em>sc</em>.
 */
// compute C = A*B (if init0==true or sc==NULL) or C += A*B (if init0==false and sc!=NULL), where
//    A is the nr-by-nrc input general matrix stored in sa[],
//    B is the nrc-by-nrc input symmetric matrix stored in packed form in sb[], and
//    C is the nr-by-nrc output general matrices stored in sc[]
// if sc==NULL is given, then the required memory will be allocated for the output matrix C;
//    otherwise it is assumed that sufficient memory has been allocated, with the address pointed to by sc
void matrix_symmatrix_multiplication(mkIndex nr, mkIndex nrc, const Real *sa, const Real *sb, Real *&sc, bool init0=true);

//! Compute <em>C=A*B</em> (if <em>init0</em>==<b>true</b> or <em>sc</em>==<tt>NULL</tt>) or <em>C+=A*B</em> (if <em>init0</em>==<b>false</b> and <em>sc</em>\f$\neq\f$<tt>NULL</tt>), where <em>A</em> is a symmetric matrix and <em>B,C</em> are general matrices.
/*! <ul>
 *  <li>  <em>A</em> is the <em>nrc</em>-by-<em>nrc</em> input symmetric matrix stored in packed form in <em>sa</em>[].
 *  <li>  <em>B</em> is the <em>nrc</em>-by-<em>nc</em> input general matrix stored in <em>sb</em>[].
 *  <li>  <em>C</em> is the <em>nrc</em>-by-<em>nc</em> output general matrix stored in <em>sc</em>[].
 *  </ul>
 *
 *  \param  nrc       (as input) is the number of rows/columns of <em>A</em>, which is also the number of rows of <em>B,C</em>.
 *  \param  nc        (as input) is the number of columns of <em>B,C</em>.
 *  \param  sa[]      (as input) stores the elements of <em>A</em>.
 *  \param  sb[]      (as input) stores the elements of <em>B</em>.
 *  \param  sc[]      (as output) stores the elements of <em>C</em>.
 *  \param  init0     (as input) tells whether to initialize <em>C</em> with zeros or not.
 *                    If <em>sc</em>==<tt>NULL</tt> is given, then <em>C</em> will be initialized with zeros anyway.
 *
 *  \remark  If <em>sc</em>==<tt>NULL</tt> is given, then the required memory of size <em>nrc*nc</em> will be allocated for the output matrix <em>C</em>;
 *           otherwise, it is assumed that sufficient memory has been allocated, with the address pointed to by <em>sc</em>.
 */
// compute C = A*B (if init0==true or sc==NULL) or C += A*B (if init0==false and sc!=NULL), where
//    A is the nrc-by-nrc input symmetric matrix stored in packed form in sa[],
//    B is the nrc-by-nc input general matrix stored in sb[], and
//    C is the nrc-by-nc output general matrix stored in sb[]
// if sc==NULL is given, then the required memory of size nrc*nc will be allocated for the output matrix C;
//    otherwise it is assumed that sufficient memory has been allocated, with the address pointed to by sc
void symmatrix_matrix_multiplication(mkIndex nrc, mkIndex nc, const Real *sa, const Real *sb, Real *&sc, bool init0=true);

/** @} */  // end of group_matrix_symmatrix



////////////////////////////////////////////////////////////////////////////
//    matrix - sparse matrix operations
////////////////////////////////////////////////////////////////////////////

/** @defgroup group_matrix_spmatrix  Matrix - Sparse Matrix Operations
 *
 *  This module contains routines `without' classes related to matrix - sparse matrix operations
 *  (declared in matkitfunc.h).
 *
 *  @{
 */

//! Compute <em>C=A+bet*B</em> (if <em>sa\f$\neq\f$sc</em>) or <em>A+=bet*B</em> (if <em>sa==sc</em>), where <em>A,C</em> are general matrices and <em>B</em> is a sparse matrix.
/*! <ul>
 *  <li>  <em>A</em> is the <em>nr</em>-by-<em>nc</em> input general matrix stored in <em>sa</em>[].
 *  <li>  <em>B</em> is the <em>nr</em>-by-<em>nc</em> input sparse matrix in CSC format stored in <em>sb</em>[], <em>rowIdx</em>[], <em>colPtr</em>[].
 *  <li>  <em>C</em> is the <em>nr</em>-by-<em>nc</em> output general matrix stored in <em>sc</em>[].
 *  </ul>
 *
 *  \param  nr        (as input) is the number of rows of <em>A,B,C</em>.
 *  \param  nc        (as input) is the number of columns of <em>A,B,C</em>.
 *  \param  sa[]      (as input) stores the elements of <em>A</em>.
 *  \param  bet       (as input) is a scalar.
 *  \param  sb[]      (as input) stores the elements of <em>B</em>.
 *  \param  rowIdx[]  (as input) stores the row indices of <em>B</em>.
 *  \param  colPtr[]  (as input) stores the column pointers of <em>B</em>.
 *  \param  sc[]      (as output) stores the elements of <em>C</em>.
 *
 *  \remark  If <em>sc</em>==<tt>NULL</tt> is given, then the required memory of size <em>nr*nc</em> will be allocated for the output matrix <em>C</em>;
 *           otherwise, it is assumed that sufficient memory has been allocated, with the address pointed to by <em>sc</em>.
 */
// compute C = A + bet*B (if sa!=sc) or A += bet*B (if sa==sc), where
//    A is the nr-by-nc input general matrix stored in sa[],
//    B is the nr-by-nc input sparse matrix in CSC format stored in sb[], rowIdx[], colPtr[], and
//    C is the nr-by-nc output general matrices stored in sc[]
// if sc==NULL is given, then the required memory of size nr*nc will be allocated for the output matrix C;
// otherwise it is assumed that sufficient memory has been allocated, with the address pointed to by sc
void matrix_spmatrix_addition(mkIndex nr, mkIndex nc, const Real *sa, Real bet, const Real *sb, const mkIndex *rowIdx, const mkIndex *colPtr, Real *&sc);

//! Compute <em>C=A*B</em> (if <em>init0</em>==<b>true</b> or <em>sc</em>==<tt>NULL</tt>) or <em>C+=A*B</em> (if <em>init0</em>==<b>false</b> and <em>sc</em>\f$\neq\f$<tt>NULL</tt>), where <em>A,C</em> are general matrices and <em>B</em> is a sparse matrix.
/*! <ul>
 *  <li>  <em>A</em> is the <em>nr</em>-by-<em>nrc</em> input general matrix stored in <em>sa</em>[].
 *  <li>  <em>B</em> is the <em>nrc</em>-by-<em>nc</em> input sparse matrix in CSC format stored in <em>sb</em>[], <em>rowIdx</em>[], <em>colPtr</em>[].
 *  <li>  <em>C</em> is the <em>nr</em>-by-<em>nc</em> output general matrix stored in <em>sc</em>[].
 *  </ul>
 *
 *  \param  nr        (as input) is the number of rows of <em>A</em>, which is also the number of rows of <em>C</em>.
 *  \param  nrc       (as input) is the number of columns of <em>A</em>, which is also the number of rows of <em>B</em>.
 *  \param  nc        (as input) is the number of columns of <em>B</em>, which is also the number of columns of <em>C</em>.
 *  \param  sa[]      (as input) stores the elements of <em>A</em>.
 *  \param  sb[]      (as input) stores the elements of <em>B</em>.
 *  \param  rowIdx[]  (as input) stores the row indices of <em>B</em>.
 *  \param  colPtr[]  (as input) stores the column pointers of <em>B</em>.
 *  \param  sc[]      (as output) stores the elements of <em>C</em>.
 *  \param  init0     (as input) tells whether to initialize <em>C</em> with zeros or not.
 *                    If <em>sc</em>==<tt>NULL</tt> is given, then <em>C</em> will be initialized with zeros anyway.
 *
 *  \remark  If <em>sc</em>==<tt>NULL</tt> is given, then the required memory of size <em>nr*nc</em> will be allocated for the output matrix <em>C</em>;
 *           otherwise, it is assumed that sufficient memory has been allocated, with the address pointed to by <em>sc</em>.
 */
// compute C = A*B (if init0==true or sc==NULL) or C += A*B (if init0==false and sc!=NULL), where
//    A is the nr-by-nrc general matrix stored in sa[],
//    B is the nrc-by-nc input sparse matrix in CSC format stored in sb[], rowIdx[], colPtr[], and
//    C is the nr-by-nc output general matrix stored in sc[]
// if sc==NULL is given, then the required memory of size nr*nc will be allocated for the output matrix C;
//    otherwise it is assumed that sufficient memory has been allocated, with the address pointed to by sc
void matrix_spmatrix_multiplication(mkIndex nr, mkIndex nrc, mkIndex nc, const Real *sa, const Real *sb, const mkIndex *rowIdx, const mkIndex *colPtr, Real *&sc, bool init0=true);

//! Compute <em>C=A*B</em> (if <em>init0</em>==<b>true</b> or <em>sc</em>==<tt>NULL</tt>) or <em>C+=A*B</em> (if <em>init0</em>==<b>false</b> and <em>sc</em>\f$\neq\f$<tt>NULL</tt>), where <em>B,C</em> are general matrices and <em>A</em> is a sparse matrix.
/*! <ul>
 *  <li>  <em>A</em> is the <em>nr</em>-by-<em>nrc</em> input sparse matrix in CSC format stored in <em>sa</em>[], <em>rowIdx</em>[], <em>colPtr</em>[].
 *  <li>  <em>B</em> is the <em>nrc</em>-by-<em>nc</em> input general matrix stored in <em>sb</em>[].
 *  <li>  <em>C</em> is the <em>nr</em>-by-<em>nc</em> output general matrix stored in <em>sc</em>[].
 *  </ul>
 *
 *  \param  nr        (as input) is the number of rows of <em>A</em>, which is also the number of rows of <em>C</em>.
 *  \param  nrc       (as input) is the number of columns of <em>A</em>, which is also the number of rows of <em>B</em>.
 *  \param  nc        (as input) is the number of columns of <em>B</em>, which is also the number of columns of <em>C</em>.
 *  \param  sa[]      (as input) stores the elements of <em>A</em>.
 *  \param  rowIdx[]  (as input) stores the row indices of <em>A</em>.
 *  \param  colPtr[]  (as input) stores the column pointers of <em>A</em>.
 *  \param  sb[]      (as input) stores the elements of <em>B</em>.
 *  \param  sc[]      (as input) stores the elements of <em>C</em>.
 *  \param  init0     (as output) tells whether to initialize <em>C</em> with zeros or not.
 *                    If <em>sc</em>==<tt>NULL</tt> is given, then <em>C</em> will be initialized with zeros anyway.
 *
 *  \remark  If <em>sc</em>==<tt>NULL</tt> is given, then the required memory of size <em>nr*nc</em> will be allocated for the output matrix <em>C</em>;
 *           otherwise, it is assumed that sufficient memory has been allocated, with the address pointed to by <em>sc</em>.
 */
// compute C = A*B (if init0==true or sc==NULL) or C += A*B (if init0==false and sc!=NULL), where
//    A is the nr-by-nrc input sparse matrix in CSC format stored in sa[], rowIdx[], colPtr[],
//    B is the nrc-by-nc input general matrix stored in sb[], and
//    C is the nr-by-nc output general matrix stored in sc[]
// if sc==NULL is given, then the required memory of size nr*nc will be allocated for the output matrix C;
//    otherwise, it is assumed that sufficient memory has been allocated, with the address pointed to by sc
void spmatrix_matrix_multiplication(mkIndex nr, mkIndex nrc, mkIndex nc, const Real *sa, const mkIndex *rowIdx, const mkIndex *colPtr,
                                    const Real *sb, Real *&sc, bool init0=true);

/** @} */  // end of group_matrix_spmatrix



////////////////////////////////////////////////////////////////////////////
//    symmetric matrix - symmetric matrix operations
////////////////////////////////////////////////////////////////////////////

/** @defgroup group_symmatrix_symmatrix  Symmetric Matrix - Symmetric Matrix Operations
 *
 *  This module contains routines `without' classes related to symmetric matrix - symmetric matrix operations
 *  (declared in matkitfunc.h).
 *
 *  @{
 */

//! Compute <em>C=A*B</em> (if <em>init0</em>==<b>true</b> or <em>sc</em>==<tt>NULL</tt>) or <em>C+=A*B</em> (if <em>init0</em>==<b>false</b> and <em>sc</em>\f$\neq\f$<tt>NULL</tt>), where <em>A,B</em> are symmetric matrices <em>C</em> is a general matrix.
/*! <ul>
 *  <li>  <em>A, B</em> are the <em>nrc</em>-by-<em>nrc</em> input symmetric matrices stored in packed form in <em>sa</em>[], <em>sb</em>[], respectively.
 *  <li>  <em>C</em> is the <em>nrc</em>-by-<em>nrc</em> output general matrix stored in <em>sc</em>[].
 *  </ul>
 *
 *  \param  nrc       (as input) is the number of rows/columns of <em>A,B,C</em>.
 *  \param  sa[]      (as input) stores the elements of <em>A</em>.
 *  \param  sb[]      (as input) stores the elements of <em>B</em>.
 *  \param  sc[]      (as output) stores the elements of <em>C</em>.
 *  \param  init0     (as input) tells whether to initialize <em>C</em> with zeros or not.
 *                    If <em>sc</em>==<tt>NULL</tt> is given, then <em>C</em> will be initialized with zeros anyway.
 *
 *  \remark  If <em>sc</em>==<tt>NULL</tt> is given, then the required memory of size <em>nrc*nrc</em> will be allocated for the output matrix <em>C</em>;
 *           otherwise, it is assumed that sufficient memory has been allocated, with the address pointed to by <em>sc</em>.
 */
// compute C = A*B (if init0==true or sc==NULL) or C += A*B (if init0==false and sc!=NULL), where
//    A, B are the nrc-by-nrc input symmetric matrices stored in packed form in sa[] and sb[], respectively, and
//    C is the nrc-by-nrc output general matrix stored in sc[]
// if sc==NULL is given, then the required memory of size nrc*nrc will be allocated for the output matrix C;
//    otherwise it is assumed that sufficient memory has been allocated, with the address pointed to by sc
void symmatrix_symmatrix_multiplication(mkIndex nrc, const Real *sa, const Real *sb, Real *&sc, bool init0=true);

/** @} */  // end of group_symmatrix_symmatrix



////////////////////////////////////////////////////////////////////////////
//    sparse matrix - sparse matrix operations
////////////////////////////////////////////////////////////////////////////

/** @defgroup group_spmatrix_spmatrix  Sparse Matrix - Sparse Matrix Operations
 *
 *  This module contains routines `without' classes related to sparse matrix - sparse matrix operations
 *  (declared in matkitfunc.h).
 *
 *  @{
 */


//! Return <b>true</b> if <em>A==B</em> mathematically, and otherwise return <b>false</b>.
/*! <ul>
 *  <li>  A is the <em>nr</em>-by-<em>nc</em> input sparse matrix in CSC format stored in <em>sa</em>[], <em>rowIdx0</em>[], colPtr0[].
 *  <li>  B is the <em>nr</em>-by-<em>nc</em> input sparse matrix in CSC format stored in <em>sb</em>[], <em>rowIdx1</em>[], colPtr1[].
 *  </ul>
 *
 *  \remark  The number of rows, <em>nr</em>, is not required in the computation, so that it is not passed.
 *  \remark  If <em>patternOnly</em>==<b>true</b>, then only the sparsity pattern is considered, and <em>sa</em>[], <em>sb</em>[] will not be used.
 *  \remark  This routine requires that the elements in each column be sorted with respect to the row indices for both <em>A</em> and <em>B</em>.
 */
// return true if A == B mathematically, and otherwise return false, where
//    A is the nr-by-nc input sparse matrix in CSC format stored in sa[], rowIdx0[], colPtr0[]
//    B is the nr-by-nc input sparse matrix in CSC format stored in sb[], rowIdx1[], colPtr1[]
// the number of rows, nr, is not required in the computation, so that it is not passed
// if patternOnly==true, then only the sparsity pattern is considered, and sa[], sb[] will not be used
// this routine requires that the elements in each column be sorted with respect to the row indices for both A and B;
bool is_spmatrix_equal_spmatrix(mkIndex nc,
     const Real *sa, const mkIndex *rowIdx0, const mkIndex *colPtr0,
     const Real *sb, const mkIndex *rowIdx1, const mkIndex *colPtr1, bool patternOnly=false);

//! Compute <em>C=alp*A+bet*B</em>, where <em>A,B,C</em> are sparse matrices and <em>alp,bet</em> are scalars.
/*! <ul>
 *  <li>  A is the <em>nr</em>-by-<em>nc</em> input sparse matrix in CSC format stored in <em>sa</em>[], <em>rowIdx0</em>[], colPtr0[].
 *  <li>  B is the <em>nr</em>-by-<em>nc</em> input sparse matrix in CSC format stored in <em>sb</em>[], <em>rowIdx1</em>[], colPtr1[].
 *  <li>  C is the <em>nr</em>-by-<em>nc</em> output sparse matrix in CSC format stored in <em>sc</em>[], <em>rowIdx2</em>[], colPtr2[].
 *  </ul>
 *
 *  \param  nc         (as input) is the number of columns of <em>A, B, C</em>.
 *  \param  alp        (as input) is a scalar.
 *  \param  sa[]       (as input) stores the elements of <em>A</em>.
 *  \param  rowIdx0[]  (as input) stores the row indices of <em>A</em>.
 *  \param  colPtr0[]  (as input) stores the column pointers of <em>A</em>.
 *  \param  bet        (as input) is a scalar.
 *  \param  sb[]       (as input) stores the elements of <em>B</em>.
 *  \param  rowIdx1[]  (as input) stores the row indices of <em>B</em>.
 *  \param  colPtr1[]  (as input) stores the column pointers of <em>B</em>.
 *  \param  sc[]       (as output) stores the elements of <em>C</em>.
 *  \param  rowIdx2[]  (as output) stores the row indices of <em>C</em>.
 *  \param  colPtr2[]  (as output) stores the column pointers of <em>C</em>.
 *  \param  patternOnly (as input) tells whether only the sparsity pattern is considered.
 *
 *  \remark  The number of rows, <em>nr</em>, is not required in the computation, so that it is not passed.
 *  \remark  If <em>sc</em>==<tt>NULL</tt> (<em>rowIdx2</em>==<tt>NULL</tt>, or <em>colPtr2</em>==<tt>NULL</tt>) is given, then sufficient memory will be allocated;
 *           otherwise, it is assumed that sufficient memory has been allocated, with the address pointed to by <em>sc</em> (<em>rowIdx2</em>, or <em>colPtr2</em>, respectively).
 *  \remark  If <em>patternOnly</em>==<b>true</b>, then only the sparsity pattern is considered, and <em>alp, bet, sa</em>[], <em>sb</em>[], and <em>sc</em>[] will not be used.
 *  \remark  This routine requires that the elements in each column be sorted with respect to the row indices for both <em>A</em> and <em>B</em>.
 *           The output matrix <em>C</em> will also have elements in each column being sorted.
 */
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
     Real *&sc, mkIndex *&rowIdx2, mkIndex *&colPtr2, bool patternOnly=false);

/** @} */  // end of group_spmatrix_spmatrix


#ifdef USE_NAMESPACE
}  // end of namespace MATKIT
#endif

#endif  // end of #ifndef MATKIT_FUNC_H
