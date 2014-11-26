//============================================================================
// the MATKIT class Vector
// coded by H.-r. Fang
// last update November, 2011
//============================================================================

#include <stdlib.h>  // for exit, atof, atoi, atol, srand, rand, etc.
#include <string.h>  // for memcpy, strcmp, strncmp, etc.
#include <stddef.h>  // for NULL pointer, pointer subtraction, etc.
#include <math.h>    // for sqrt, pow, etc.
#include <iostream>  // for cout, cerr, endl, etc. (under namespace std)

#include "matkit.h"

using std::endl;

#ifdef USE_NAMESPACE
namespace MATKIT {
#endif


// display vector v with operator <<
// in other words, print out vector v to ostream s
// this routine invokes print_vector()
std::ostream& operator<<(std::ostream &s, const Vector &v) {
    print_vector(s, v.Length(), v.Store());
    return s;
}

// read a vector from a plain text file fileName[] and return the result
// this member function invokes vector_read(); see the comments therein for the file format
Vector Vector::read(const char fileName[]) {
    mkIndex len;
    Real *sv = NULL;

    // read the vector
    int info = vector_read(fileName, len, sv);

    // in case of any error, report and quit
    switch (info) {
        case -1:  // file open error
            *Basics::err << "Vector::read(const char*, bool): file open error!" << endl;
            Basics::quit(1);
        case -2:  // file is empty
            *Basics::err << "Vector::read(const char*): file is empty!" << endl;
            Basics::quit(1);
        case -3:  // length of the vector error
            *Basics::err << "Vector::read(const char*): length of the vector error!" << endl;
            Basics::quit(1);
        case -4:  // length of the vector error
            *Basics::err << "Vector::read(const char*): the vector data is incomplete or invalid!" << endl;
            Basics::quit(1);
        default:  // a successful read
            if (info != 0) {
                *Basics::err << "Vector::read(const char*): invalid return flag of vector_read()!" << endl;
                Basics::quit(1);
            }
    }

    return Vector(len, sv);
}


// return a vector of zeros of length n
// this is a static function
Vector Vector::zeros(mkIndex n) {
    // allocate memory
    Real *store = new Real[n];
    // memset(store, 0, n*sizeof(Real));  // works for IEEE 754, not guaranteed working for other types of floating point numbers

    // set it zero
    Real *sv = store;
    mkIndex ii = n;
    while (ii--)
        *sv++ = 0.0;

    // return the result
    return Vector(n, store);
}

// return a vector of ones of length n
// this is a static function
Vector Vector::ones(mkIndex n) {
    // allocate memory
    Real *store = new Real[n];

    // set ones
    Real *sv = store;
    mkIndex ii = n;
    while (ii--)
        *sv++ = 1.0;

    // return the result
    return Vector(n, store);
}

// return a vector of length n containing pseudo-random values drawn from a uniform distribution on the interval [a,b]
// this is a static function
Vector Vector::random(mkIndex n, Real a, Real b) {
    // exception handling
    if (n <= 0)
        return Vector();

    // allocate memory
    Real *store = new Real[n];

    // set random elements
    Real *sv = store;
    mkIndex ii = n;
    while (ii--)
        *sv++ = a+(b-a)*((Real)rand()/(Real)RAND_MAX);

    // return the result
    return Vector(n, store);
}

// concatenate two vectors v and w
// this static function invokes vector_concatenate()
Vector Vector::concatenate(const Vector &v, const Vector &w) {
    // pointer to the storage of the output vector
    Real *su = NULL;

    // catenate v and w
    vector_concatenate(v.Length(), v.Store(), w.Length(), w.Store(), su);

    // return the result
    return Vector(v.Length()+w.Length(), su);
}

// inner product v'*w
// this static function invokes vector_inner_product()
Real Vector::innerProduct(const Vector &v, const Vector &w) {
    // compute v(1)*w(1) + v(2)*w(2) + ... + v(len)*w(len)

    // exception handling
    mkIndex len = v.Length();
    if (len != w.Length()) {
        *Basics::err << "InnerProduct(const Vector&, const Vector&): incompatible dimensions!" << endl;
        Basics::quit(1);
    }
    // compute the inner product
#ifdef USE_MWBLAS
    // invoke xDOT
    // MWBLAS routines sdot and ddot require non-const storage for the input vectors,
    // so it cannot be incorporated into vector_inner_product and therefore is put here
    mwSignedIndex len2 = (mwSignedIndex)len;
    mwSignedIndex inc = 1;
    #ifdef USE_SINGLE
        return sdot(&len2, v.store, &inc, w.store, &inc);
    #else
        return ddot(&len2, v.store, &inc, w.store, &inc);
    #endif
#else
    return vector_inner_product(len, v.Store(), w.Store());
#endif
}

// (scaled) outer product alp*v*w'
// this static function invokes vector_outer_product()
Matrix Vector::outerProduct(const Vector &v, const Vector &w, Real alp) {
   // get the dimensions
    mkIndex nrows = v.Length();
    mkIndex ncols = w.Length();

#ifdef USE_MWBLAS
    // MWBLAS routines sger and dger require non-const storage for the input vectors,
    // so it cannot be incorporated into vector_outer_product and therefore is put here

    // allocate memory
    mkIndex sz = nrows*ncols;
    Real *sa = new Real[sz];

    // zero initialization
    Real *a0 = sa;
    while (sz--)
        *a0++ = 0.0;

    // invoke xGER
    mwSignedIndex nrows2 = (mwSignedIndex)nrows;
    mwSignedIndex ncols2 = (mwSignedIndex)ncols;
    mwSignedIndex inc = 1;
    #ifdef USE_SINGLE
        sger(&nrows2, &ncols2, &alp, v.store, &inc, w.store, &inc, sa, &nrows2);
    #else
        dger(&nrows2, &ncols2, &alp, v.store, &inc, w.store, &inc, sa, &nrows2);
    #endif
#else
    Real *sa = NULL;
    vector_outer_product(nrows, ncols, alp, v.Store(), w.Store(), sa);
#endif

    // return the result
    return Matrix(nrows, ncols, sa);
}

// (scaled) outer product alp*v*v'
// this static function invokes vector_outer_product()
SymmetricMatrix Vector::outerProduct(const Vector &v, Real alp) {
    // get the sizes
    mkIndex nrc = v.Length();

    // compute the outer product alp*v*v'
#ifdef USE_MWBLAS
    // MWBLAS routines sspr and dspr require non-const storage for the input vector,
    // so it cannot be incorporated into vector_outer_product and therefore is put here

    // allocate memory
    mkIndex sz = nrc*(nrc+1u)/2u;
    Real *sa = new Real[sz];

    // zero initialization
    Real *a0 = sa;
    while (sz--)
        *a0++ = 0.0;

    // invoke xSPR
    char uplo = 'L';
    mwSignedIndex inc = 1;
    mwSignedIndex nrc2 = (mwSignedIndex)nrc;
    #ifdef USE_SINGLE
        sspr(&uplo, &nrc2, &alp, v.store, &inc, sa);
    #else
        dspr(&uplo, &nrc2, &alp, v.store, &inc, sa);
    #endif
#else
    Real *sa = NULL;
    vector_outer_product(nrc, alp, v.Store(), sa);
#endif

    // return the result
    return SymmetricMatrix(nrc, sa);
}

// (scaled) symmetric outer product alp*(v*w'+w*v')
// this static function invokes vector_symmetric_outer_product()
SymmetricMatrix Vector::symmetricOuterProduct(const Vector &v, const Vector &w, Real alp) {
    // check the dimension
    mkIndex nrc = v.Length();
    if (nrc != w.Length()) {
        *Basics::err << "SymmetricOuterProduct(const Vector&, const Vector&, Real): vector dimensions must agree!" << endl;
        Basics::quit(1);
    }

    // compute the symmetric outer product
#ifdef USE_MWBLAS
    // allocate memory
    mkIndex sz = nrc*(nrc+1u)/2u;
    Real *sa = new Real[sz];

    // zero initialization
    Real *a0 = sa;
    while (sz--)
        *a0++ = 0.0;

    // invoke xSPR2
    char uplo = 'L';
    mwSignedIndex inc = 1;
    mwSignedIndex nrc2 = (mwSignedIndex)nrc;
    #ifdef USE_SINGLE
        sspr2(&uplo, &nrc2, &alp, v.store, &inc, w.store, &inc, sa);
    #else
        dspr2(&uplo, &nrc2, &alp, v.store, &inc, w.store, &inc, sa);
    #endif
#else
    Real *sa = NULL;
    vector_symmetric_outer_product(nrc, alp, v.Store(), w.Store(), sa);
#endif

    // return the result
    return SymmetricMatrix(nrc, sa);
}


// constructors
// a constructor for an empty vector
Vector::Vector() {
    // empty vector
    length = 0;
    store = NULL;
}

// a constructor for a vector of length n of zeros
Vector::Vector(mkIndex n) {
    // set the length
    length = n;

    // allocate memory
    if (n > 0)
        store = new Real[n];
    else
        store = NULL;

    // zero initialization
    if (n > 0) {
        Real *v0 = store;
        while (n--)
            *v0++ = 0.0;
    }
}

// a constructor for a vector v of length n with elements stored in store0[] in order
// note that it performs a shallow copy of store0[] and the memory will be freed by the destructor (i.e. `delete [] store0')
Vector::Vector(mkIndex n, Real *store0) {
    // set the length
    length = n;

    // a shallow copy
    store = store0;
}

// a copy constructor
// this constructor invokes memory_xcopy()
Vector::Vector(const Vector &v) {
    // set length
    length = v.Length();

    // allocate memory
    store = new Real[length];

    // copy elements
    memory_xcopy(v.Length(), v.Store(), store);
}

// a destructor
Vector::~Vector() {
    delete [] store;
}


// resize *this to be of length n
void Vector::resize(mkIndex n) {
    if (n != length) {
        // free memory and allocate memory precisely required
        delete [] store;
        if (n > 0)
            store = new Real[n];
        else
            store = NULL;
        length = n;
    }
}

// set *this be a vector of length n and stored in store0[]
// it performs a shallow copy of store0[] and the memory will be freed by the destructor (i.e. `delete [] store0')
void Vector::set(mkIndex n, Real *store0) {
    length = n;
    if (store != store0) {
        // free memory and do a shallow copy
        delete [] store;
        store = store0;
    }
}

// return a (square) matrix of dimension len+|k| with the k-th diagonal formed by the elements of *this in order, and otherwise zero
// this member function invokes vector_to_diag()
// this member function resembles diag(v) in OCTAVE/MATLAB
Matrix Vector::diag(mkSignedIndex k) const {
    // dimension of the output matrix
    mkIndex nrc = (k>=0) ? (length+k) : (length-k);

    // pointer for output
    Real *sa = NULL;

    // form the matrix
    vector_to_diag(length, Store(), sa, k);

    // return the result
    return Matrix(nrc, nrc, sa);
}

// return a (square) sparse matrix of dimension len+|k| with the k-th diagonal formed by the elements of *this in order, and otherwise zero
// this member function invokes vector_to_spdiag()
// this member function resembles spdiag(v) in OCTAVE/MATLAB
SparseMatrix Vector::spdiag(mkSignedIndex k) const {
    // dimension of the output matrix
    mkIndex nrc = (k>=0) ? (length+k) : (length-k);

    // pointers for the output matrix
    Real *sa = NULL;
    mkIndex *rowIdx = NULL;
    mkIndex *colPtr = NULL;

    // form the matrix
    vector_to_spdiag(length, store, sa, rowIdx, colPtr, k);

    // return the result
    return SparseMatrix(nrc, nrc, sa, rowIdx, colPtr);
}


// a subvector formed by elements i,...,j of *this
// this member function invokes memory_xcopy()
Vector Vector::subVector(mkIndex i, mkIndex j) const {
    if (i == 0 || j == 0) {
        *Basics::err << "Vector::subVector(mkIndex, mkIndex) const: subscript indices must be positive!" << endl;
        Basics::quit(1);
    }
    if (i > length || j > length) {
        *Basics::err << "Vector::subVector(mkIndex, mkIndex) const: index exceeds vector length!" << endl;
        Basics::quit(1);
    }

    // length of the output vector
    mkIndex len = (i<=j) ? (j-i+1u) : (i-j+1u);

    // allocate memory
    Real *sw = NULL;
    if (len)
        sw = new Real[len];

    // copy v(i:j) if i<=j, or v(i:-1:j) if i>j, where v is *this
    mkSignedIndex incx = (i<=j) ? 1 : -1;
    const Real *a0 = (i<=j) ? Store()+(i-1u) : Store()+(j-1u);
    memory_xcopy(len, a0, sw, incx, 1);

    // return the result
    return Vector(len, sw);
}


// 2-norm of *this
Real Vector::norm2() const {
    // compute square root of v(1)^2 + v(2)^2 + ... + v(length)^2, where v is *this
#if defined(USE_MWBLAS)
    mwSignedIndex length2 = (mwSignedIndex)length;
    mwSignedIndex inc = 1;
    #ifdef USE_SINGLE
        return snrm2(&length2, store, &inc);
    #else
        return dnrm2(&length2, store, &inc);
    #endif
#elif defined(USE_BLAS)
    int length2 = (int)length;
    int inc = 1;
    #ifdef USE_SINGLE
        return snrm2_(&length2, store, &inc);
    #else
        return dnrm2_(&length2, store, &inc);
    #endif
#else
    return sqrt(sumSquare());
#endif
}

// Frobenius norm of *this, the same as 2-norm since it is a vector
Real Vector::normFrobenius() const {
    return norm2();
}

// sum of squared elements
Real Vector::sumSquare() const {
    // compute v(1)^2 + v(2)^2 +...+ v(length)^2, where v is *this
    Real sum = 0.0;
    Real *s = store;
    mkIndex ii = length;
    while (ii--)
        sum += Basics::square(*s++);

    // return the result
    return sum;
}

// 1-norm of *this
Real Vector::norm1() const {
    // compute t
    Real nrm = 0.0;
    Real *ss = store;
    mkIndex ii = length;
    while (ii--)
        nrm += Basics::abs(*ss++);

    // return the result
    return nrm;
}

// infinity-norm of *this
Real Vector::normInfinity() const {
    // compute the norm
    Real *ss = store;
    Real nrm = 0.0;
    mkIndex ii = length;
    while (ii--) {
        nrm = (nrm>=(*ss) ? nrm : (*ss));
        ss ++;
    }

    // return the result
    return nrm;
}


// access v(i)
// return a reference to v(i), i.e. can be used to write v(i), where v is *this
Real& Vector::operator()(mkIndex i) {
    // exception handling
    if (i == 0) {
        *Basics::err << "Vector::operator()(mkIndex): subscript indices must be positive!" << endl;
        Basics::quit(1);
    }
    if (i > length) {
        *Basics::err << "Vector::operator()(mkIndex): index exceeds vector length!" << endl;
        Basics::quit(1);
    }

    // return a reference to v(i)
    return store[i-1u];
}

// read v(i) via a const reference or pointer to v, where v is *this
Real Vector::operator()(mkIndex i) const {
    // exception handling
    if (i == 0) {
        *Basics::err << "Vector::operator()(mkIndex) const: subscript indices must be positive!" << endl;
        Basics::quit(1);
    }
    if (i > length) {
        *Basics::err << "Vector::operator()(mkIndex) const: index exceeds vector length!" << endl;
        Basics::quit(1);
    }

    // read v(i), via a const reference or pointer
    return store[i-1u];
}


// operator =
// return a const reference to v after setting v(i)=f for i=1,...,length, where v is *this
const Vector& Vector::operator=(Real f) {
    // set elements be f
    mkIndex ii = length;
    Real *s = store;
    while (ii--)
        *s++ = f;

    // return the result
    return *this;
}

// return a const reference to w after setting w=v, where w is *this
// this operator invokes memory_xcopy()
const Vector& Vector::operator=(const Vector &v) {
    // length of the vector
    mkIndex len = v.Length();

    if (length != len) {
        // free memory and allocated memory precisely required
        delete [] store;
        if (len)
            store = new Real[len];
        else
            store = NULL;
        length = len;
    }

    // copy elements
    memory_xcopy(len, v.Store(), store);

    // return the result
    return *this;
}

// operator -
// return -w, where w is *this
// this operator invokes vector_scalar_operation()
Vector Vector::operator-() const {
    // pointer for output
    Real *sv = NULL;

    // elementwise multiply -1.0
    vector_scalar_operation(length, Store(), -1.0, sv, 2);

    // return the result
    return Vector(length, sv);
}


// operators +,-,*,/,+=,-=,*=,/= a real number
// return v with v(i)=w(i)+f for i=1,...,length, where w is *this
// this operator invokes vector_scalar_operation()
Vector Vector::operator+(Real f) const {
    // pointer for output
    Real *sv = NULL;

    // compute v = *this + f, with v stored in sv[]
    vector_scalar_operation(length, Store(), f, sv, 0);

    // return the result
    return Vector(length, sv);
}

// return v with v(i)=w(i)-f for i=1,...,length, where w is *this
// this operator invokes vector_scalar_operation()
Vector Vector::operator-(Real f) const {
    // pointer for output
    Real *sv = NULL;

    // compute v = *this + f, with v stored in sv[]
    vector_scalar_operation(length, Store(), f, sv, 1);

    // return the result
    return Vector(length, sv);
}

// return v with v(i)=w(i)*f for i=1,...,length, where w is *this
// this operator invokes vector_scalar_operation()
Vector Vector::operator*(Real f) const {
    // pointer for output
    Real *sv = NULL;

    // compute v = *this + f, with v stored in sv[]
    vector_scalar_operation(length, Store(), f, sv, 2);

    // return the result
    return Vector(length, sv);
}

// return v with v(i)=w(i)/f for i=1,...,length, where w is *this
// this operator invokes vector_scalar_operation()
Vector Vector::operator/(Real f) const {
    // exception handling
    if (f == 0.0) {
        *Basics::err << "Vector::operator/(Real) const: division-by-zero exception!" << endl;
        Basics::quit(1);
    }
    // pointer for output
    Real *sv = NULL;

    // compute v = *this + f, with v stored in sv[]
    vector_scalar_operation(length, Store(), f, sv, 3);

    // return the result
    return Vector(length, sv);
}

// return a const reference to w after computing w(i)+=f for i=1,...,length, where w is *this
// this operator invokes vector_scalar_operation()
const Vector& Vector::operator+=(Real f) {
    // compute *this += f
    vector_scalar_operation(length, Store(), f, store, 0);

    // return the result
    return *this;
}

// return a const reference to w after computing w(i)-=f for i=1,...,length, where w is *this
// this operator invokes vector_scalar_operation()
const Vector& Vector::operator-=(Real f) {
    // compute *this -= f
    vector_scalar_operation(length, Store(), f, store, 1);

    // return the result
    return *this;
}

// return a const reference to w after computing w(i)*=f for i=1,...,length, where w is *this
// this operator invokes vector_scalar_operation()
const Vector& Vector::operator*=(Real f) {
    // compute v *= f (elementwise), where v is *this
    vector_scalar_operation(length, Store(), f, store, 2);

    // return the result
    return *this;
}

// return a const reference to w after computing w(i)/=f for i=1,...,length, where w is *this
// this operator invokes vector_scalar_operation()
const Vector& Vector::operator/=(Real f) {
    // exception handling
    if (f == 0.0) {
        *Basics::err << "Vector::operator/=(Real): division-by-zero exception!" << endl;
        Basics::quit(1);
    }
    // compute *this /= f
    vector_scalar_operation(length, Store(), f, store, 3);

    // return the result
    return *this;
}


// a real number +,-,* a vector
// return w with w(i)=f+v(i) for i=1,...,v.length
// this operator (eventually) invokes vector_scalar_operation()
Vector operator+(Real f, const Vector &v) {
    return v+f;
}

// return w with w(i)=f-v(i) for i=1,...,v.length
// this operator (eventually) invokes vector_scalar_operation()
Vector operator-(Real f, const Vector &v) {
    return -v+f;
}

// return f*v
// this operator (eventually) invokes vector_scalar_operation()
Vector operator*(Real f, const Vector &v) {
    return v*f;
}


// vector-vector addition/subtraction
// return w + v, where w is *this
// this operator invokes elementwise_addition()
Vector Vector::operator+(const Vector &v) const {
    // exception handling
    if (length != v.Length()) {
        *Basics::err << "Vector::operator+(const Vector&) const: incompatible dimensions!" << endl;
        Basics::quit(1);
    }

    // pointer for output
    Real *su = NULL;

    // compute u = *this + v, with u stored in su[]
    elementwise_addition(length, Store(), 1.0, v.Store(), su);

    // return the result
    return Vector(length, su);
}

// return w - v, where w is *this
// this operator invokes elementwise_addition()
Vector Vector::operator-(const Vector &v) const {
    // exception handling
    if (length != v.Length()) {
        *Basics::err << "Vector::operator-(const Vector&) const: incompatible dimensions!" << endl;
        Basics::quit(1);
    }

    // pointer for output
    Real *su = NULL;

    // compute u = *this - v, with u stored in su[]
    elementwise_addition(length, Store(), -1.0, v.Store(), su);

    // return the result
    return Vector(length, su);
}

// return a const reference to w after computing w = w+v, where w is *this
// this operator invokes elementwise_addition()
const Vector& Vector::operator+=(const Vector &v) {
    // exception handling
    if (length != v.Length()) {
        *Basics::err << "Vector::operator+=(const Vector&): incompatible dimensions!" << endl;
        Basics::quit(1);
    }

    // compute *this += v
    elementwise_addition(length, Store(), 1.0, v.Store(), store);

    // return the result
    return *this;
}

// return a const reference to w after computing w = w-v, where w is *this
// this operator invokes elementwise_addition()
const Vector& Vector::operator-=(const Vector &v) {
    // exception handling
    if (length != v.Length()) {
        *Basics::err << "Vector::operator-=(const Vector&): incompatible dimensions!" << endl;
        Basics::quit(1);
    }

    // compute *this -= v
    elementwise_addition(length, Store(), -1.0, v.Store(), store);

    // return the result
    return *this;
}


// vector - matrix multiplication
// compute v' = w'*A and return the result v, where w is *this
// this operator invokes vector_matrix_multiplication()
Vector Vector::operator*(const Matrix &A) const {
    // exception handling
    mkIndex nr = A.Nrows();
    mkIndex nc = A.Ncols();
    if (length != nr) {
        *Basics::err << "Vector Vector::operator*(const Matrix &) const: inner vector - matrix dimensions must agree!" << endl;
        Basics::quit(1);
    }
#ifdef USE_MWBLAS
    // MWBLAS routines sgemv and dgemv require non-const storage for the input vector and matrix,
    // so it cannot be incorporated into matrix_vector_multiplication and therefore is put here
    Real *sv = new Real[nc];
    char trans = 'T';
    Real alp = 1.0, bet = 0.0;
    mwSignedIndex inc = 1;
    mwSignedIndex nrows2 = (mwSignedIndex)nr;
    mwSignedIndex ncols2 = (mwSignedIndex)nc;
    #ifdef USE_SINGLE
        sgemv(&trans, &nrows2, &ncols2, &alp, A.store, &nrows2, store, &inc, &bet, sv, &inc);
    #else
        dgemv(&trans, &nrows2, &ncols2, &alp, A.store, &nrows2, store, &inc, &bet, sv, &inc);
    #endif
#else
    Real *sv = NULL;
    vector_matrix_multiplication(A.Nrows(), A.Ncols(), Store(), A.Store(), sv);
#endif

    // return the result
    return Vector(nc, sv);
}

// compute w' = w'*A and return the const reference to w, where w is *this
// this operator (eventually) invokes vector_matrix_multiplication()
const Vector& Vector::operator*=(const Matrix &A) {
    // exception handling
    if (A.Nrows() != length) {
        *Basics::err << "Vector::operator*=(const Matrix&): inner vector - matrix dimensions must agree!" << endl;
        Basics::quit(1);
    }

    // self-explained
    return *this = (*this)*A;
}


// vector - symmetric matrix multiplication
// compute v' = w'*A and return the result v, where w is *this
// this operator (eventually) invokes symmatrix_vector_multiplication()
Vector Vector::operator*(const SymmetricMatrix &A) const {
    // exception handling
    if (A.Nrows() != length) {
        *Basics::err << "Vector::operator*(const SymmetricMatrix&): inner vector - matrix dimensions must agree!" << endl;
        Basics::quit(1);
    }

    // self-explained; not that A is symmetric
    return A*(*this);
}

// compute w' = w'*A and return the const reference to w, where w is *this
// this operator (eventually) invokes symmatrix_vector_multiplication()
const Vector& Vector::operator*=(const SymmetricMatrix &A) {
    // exception handling
    if (A.Nrows() != length) {
        *Basics::err << "Vector::operator*=(const SymmetricMatrix&): inner vector - matrix dimensions must agree!" << endl;
        Basics::quit(1);
    }

    // self-explained; note that A is symmetric
    return *this = A*(*this);
}


// vector - sparse matrix multiplication
// compute v' = w'*A and return the result v, where w is *this
// this operator invokes vector_spmatrix_multiplication()
Vector Vector::operator*(const SparseMatrix &A) const {
    // exception handling
    if (A.Nrows() != length) {
        *Basics::err << "Vector::operator*(const SparseMatrix&): inner vector - matrix dimensions must agree!" << endl;
        Basics::quit(1);
    }
    // pointer for output
    Real *sv = NULL;

    // compute v' = w'*A, where w is *this and v is stored in sv[]
    vector_spmatrix_multiplication(A.Ncols(), store, A.Store(), A.RowIndex(), A.ColPointer(), sv);  // declared in spmatrix.h and implemented in spmatrix.cpp

    // return the result
    return Vector(A.Ncols(), sv);
}

// compute w' = w'*A and return the const reference to w, where w is *this
// this operator (eventually) invokes vector_spmatrix_multiplication()
const Vector& Vector::operator*=(const SparseMatrix &A) {
    // exception handling
    if (A.Nrows() != length) {
        *Basics::err << "Vector::operator*=(const SparseMatrix&): inner vector - matrix dimensions must agree!" << endl;
        Basics::quit(1);
    }

    // self-explained
    return *this = (*this)*A;
}


// return a vector with elements as those of *this reversed
// to be precise, return v such that v(i) = w(length-i+1) for i=1,...,length, where w is *this
Vector Vector::reverse() const {
    // allocate memory
    Real *sv = new Real[length];

    // copy elements reversely
    Real *v0 = sv+length;
    const Real *w0 = store;
    mkIndex ii = length;
    while (ii--)
        *(--v0) = *(w0++);

    // return the result
    return Vector(length, sv);
}


// orthogonalize v against Q(:,1),...,Q(:,ncols), where
// v is *this and ncols is the number of columns in Q, assuming
// Q(:,i) and Q(:,j) are orthogonal to each other for i != j
// if normalize==true, then v will be normalized (i.e. v = v/||v||)
void Vector::modifiedGramSchmidt(const Matrix &Q, bool normalize) {
    // exception handling
    if (Q.Nrows() != length) {
        *Basics::err << "Vector::modifiedGramSchmit(const Matrix&, bool): incompatible dimensions!" << endl;
        Basics::quit(1);
    }

    const Real *sq = Q.Store();
    mkIndex nc = Q.Ncols();
    while (nc--) {
        // orthogonalize v against the i-th column of Q, where
        // v is *this and i = Q.Ncols()-nc
        // in the words, compute v = v - <v,q>*q, where
        // q is Q(:,i) and <v,q> is the inner product of v and q

        // sq currently points to Q(1,i)

        // compute inner_product = <v,q>
        Real inner_product = vector_inner_product(length, Store(), sq);

        // compute v = v - inner_product * q
        elementwise_addition(length, Store(), -inner_product, sq, store);

        // for pointer to Q(1,i+1)
        sq += Q.Nrows();
    }

    if (normalize) {
        // normalize the result
        *this /= this->norm2();
    }
}

// orthogonalize v against Q(:,k) for flag[k-1]==true, k=1,...,ncols, where
// v is *this and len is the number of columns in Q
// it is assumed that Q(:,i) and Q(:,j) are orthogonal to each other for i != j
// if normalize==true, then v will be normalized (i.e. v = v/||v||)
void Vector::modifiedGramSchmidt(const Matrix &Q, const bool *flag, bool normalize) {
    // exception handling
    if (Q.Nrows() != length) {
        *Basics::err << "Vector::modifiedGramSchmit(const Matrix&, bool*, bool): incompatible dimensions!" << endl;
        Basics::quit(1);
    }

    if (flag == NULL) {
        // orthogonalize Q(:,i) for all i=1,...,ncols, with ncols the number of columns in Q
        this->modifiedGramSchmidt(Q, normalize);
        return;
    }

    const Real *sq = Q.Store();
    mkIndex nc = Q.Ncols();
    const bool *flg = flag;
    while (nc--) {
        if (*flg++) {
            // orthogonalize v against Q(:,i), where
            // v is *this and i = Q.Ncols()-nc
            // in the words, compute v = v - <v,q>*q, where
            // q is Q(:,i) and <v,q> is the inner product of v and q

            // sq currently points to Q(1,i)

            // compute inner_product = <v,q>
            Real inner_product = vector_inner_product(length, Store(), sq);

            // compute v = v - inner_product * q
            elementwise_addition(length, Store(), -inner_product, sq, store);

            // for pointer to Q(1,i+1)
            sq += Q.Nrows();
        }
        // else Q(:,i) is not flagged, so do not orthogonalize against it
    }

    if (normalize) {
        // normalize the result
        *this /= this->norm2();
    }
}

#ifdef USE_MEX
// convert *mxArray to Vector
Vector Vector::mxArray2Vector(const mxArray *mv) {
    // a mkIndex and a Real pointer for the output vector
    mkIndex len;
    Real *sv = NULL;

    // convert
    mxArray_to_vector(mv, len, sv);

    // return the result
    return Vector(len, sv);
}

// convert Vector to *mxArray
mxArray *Vector::Vector2mxArray(const Vector &v) {
    return vector_to_mxArray(v.Length(), v.Store());
}
#endif  // of #ifdef USE_MEX


#ifdef USE_NAMESPACE
}  // end of namespace NEWMAT
#endif
