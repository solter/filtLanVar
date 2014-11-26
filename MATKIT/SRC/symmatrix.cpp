//============================================================================
// the MATKIT class SymmetricMatrix
// coded by H.-r. Fang
// last update August, 2011
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


// display symmetric matrix A with operator <<
// in other words, print out symmetric matrix A to ostream s
// this routine invokes print_symmatrix()
std::ostream& operator<<(std::ostream &s, const SymmetricMatrix &A) {
    print_symmatrix(s, A.Nrows(), A.Store());
    return s;
}


// return an n-by-n (symmetric) matrix of zeros
// this is a static function
SymmetricMatrix SymmetricMatrix::zeros(mkIndex n) {
    // allocate memory
    mkIndex sz = n*(n+1u)/2u;
    Real *sa = new Real[sz];
    // memset(sa, 0, sizeof(Real)*n*(n+1)/2);  // works for IEEE 754, not guaranteed working for other types of floating point numbers

    // set it zero
    Real *a0 = sa;
    while (sz--)
        *a0++ = 0.0;

    // return the result
    return SymmetricMatrix(n, sa);
}

// return an n-by-n (symmetric) matrix of ones
// this is a static function
SymmetricMatrix SymmetricMatrix::ones(mkIndex n) {
    // allocate memory
    mkIndex sz = n*(n+1u)/2u;
    Real *sa = new Real[sz];

    // set ones
    Real *a0 = sa;
    while (sz--)
        *a0++ = 1.0;

    // return the result
    return SymmetricMatrix(n, sa);
}

// return an n-by-n (symmetric) matrix with ones on the main diagonal and zeros elsewhere
// this is a static function
SymmetricMatrix SymmetricMatrix::eye(mkIndex n) {
    // allocate memory
    mkIndex sz = n*(n+1u)/2u;
    Real *sa = new Real[sz];

    // set the eye
    Real *a0 = sa;
    mkIndex nn = n;
    while (nn--) {
        *a0 = 1.0;
        mkIndex ii = nn;
        while (ii--)
            *a0 = 0.0;
    }

    // return the result
    return SymmetricMatrix(n, sa);
}

// return an n-by-n symmetric matrix containing pseudo-random values drawn from a uniform distribution on the interval [a,b]
// this is a static function
SymmetricMatrix SymmetricMatrix::random(mkIndex n, Real a, Real b) {
    // allocate memory
    mkIndex sz = n*(n+1u)/2u;
    Real *sa = new Real[sz];

    // set random elements
    Real *a0 = sa;
    while (sz--)
        *a0++ = a+(b-a)*((Real)rand()/(Real)RAND_MAX);

    // return the result
    return SymmetricMatrix(n, sa);
}


// a constructor for an empty (symmetric) matrix
SymmetricMatrix::SymmetricMatrix() {
    nrows = 0;
    store = NULL;
}

// a constructor for an n-by-n (symmetric) matrix of zeros
SymmetricMatrix::SymmetricMatrix(mkIndex n) {
    // set dimension
    nrows = n;

    // allocate memory
    mkIndex sz = n*(n+1u)/2u;
    store = new Real[sz];

    // set it zero
    Real *sa = store;
    while (sz--)
        *sa++ = 0.0;
}

// a constructor for an n-by-n symmetric matrix of with elements stored columnwise in packed form in store0[]
// note that it performs a shallow copy of store0 and the memory will be freed by the destructor (i.e. `delete [] store0')
SymmetricMatrix::SymmetricMatrix(mkIndex n, Real *store0) {
    // set dimension
    nrows = n;

    // a shallow copy
    store = store0;
}

// a copy constructor
// this constructor invokes memory_xcopy()
SymmetricMatrix::SymmetricMatrix(const SymmetricMatrix &B) {
    // set dimension
    nrows = B.Nrows();

    // allocate memory
    mkIndex sz = nrows*(nrows+1u)/2u;
    store = new Real[sz];

    // copy elements
    memory_xcopy(sz, B.Store(), store);
}


// a destructor
SymmetricMatrix::~SymmetricMatrix() {
    delete [] store;
}


// convert *this to a general matrix
// this member function invokes symmetric_to_general()
Matrix SymmetricMatrix::toGeneral() const {
    // pointer for output
    Real *sa = NULL;

    // convert *this to a general matrix stored in sa[]
    symmetric_to_general(nrows, Store(), sa);

    // return the result
    return Matrix(nrows, nrows, sa);
}

// resize *this to be of dimensions nrc-by-nrc
void SymmetricMatrix::resize(mkIndex nrc) {
    if (nrc != nrows) {
        // free memory and allocate the memory precisely required
        delete [] store;
        store = new Real[nrc*(nrc+1u)/2u];

        // set dimension
        nrows = nrc;
    }
}

// set *this as an nrc-by-nrc symmetric matrix stored in store0[] columnwise in packed form
// it performs a shallow copy of store0[] and the memory will be freed by the destructor (i.e. `delete [] store0')
void SymmetricMatrix::set(mkIndex nrc, Real *store0) {
    // set dimension
    nrows = nrc;

    if (store != store0) {
        // free memory and do a shallow copy
        delete [] store;
        store = store0;
    }
}


// vector formed by column j (or equivalently row j)
// this member function invokes elements_of_symmatrix()
Vector SymmetricMatrix::column(mkIndex j) const {
    if (j == 0) {
        *Basics::err << "SymmetricMatrix::column(mkIndex) const: subscript indices must be positive!" << endl;
        Basics::quit(1);
    }
    if (j > nrows) {
        *Basics::err << "SymmetricMatrix::column(mkIndex) const: index exceeds matrix dimensions!" << endl;
        Basics::quit(1);
    }

    // extract the column
    Real *sv = NULL;
    elements_of_symmatrix(nrows, store, j, 1u, nrows, sv);

    // return the result
    return Vector(nrows, sv);
}

// vector formed by row i (or equivalently column i)
// this member function invokes elements_of_symmatrix()
Vector SymmetricMatrix::row(mkIndex i) const {
    if (i == 0) {
        *Basics::err << "SymmetricMatrix::row(mkIndex) const: subscript indices must be positive!" << endl;
        Basics::quit(1);
    }
    if (i > nrows) {
        *Basics::err << "SymmetricMatrix::row(mkIndex) const: index exceeds matrix dimensions!" << endl;
        Basics::quit(1);
    }

    // extract the row
    Real *sv = NULL;
    elements_of_symmatrix(nrows, store, i, 1u, nrows, sv);

    // return the result
    return Vector(nrows, sv);
}

// submatrix formed by columns j1,...,j2
// this member function invokes submatrix_of_symmatrix()
Matrix SymmetricMatrix::columns(mkIndex j1, mkIndex j2) const {
    if (j1 == 0 || j2 == 0) {
        *Basics::err << "SymmetricMatrix::columns(mkIndex, mkIndex) const: subscript indices must be positive!" << endl;
        Basics::quit(1);
    }
    if (j1 > nrows || j2 > nrows) {
        *Basics::err << "SymmetricMatrix::columns(mkIndex, mkIndex) const: index exceeds matrix dimensions!" << endl;
        Basics::quit(1);
    }

    // find dimensions
    mkIndex nc = (j1<=j2) ? (j2-j1+1u) : (j1-j2+1u);

    // output storage address
    Real *sb = NULL;

    // copy elements
    submatrix_of_symmatrix(nrows, store, 1, nrows, j1, j2, sb);

    // return the result
    return Matrix(nrows, nc, sb);
}

// submatrix formed by columns i1,...,i2
// this member function invokes submatrix_of_symmatrix()
Matrix SymmetricMatrix::rows(mkIndex i1, mkIndex i2) const {
    if (i1 == 0 || i2 == 0) {
        *Basics::err << "SymmetricMatrix::rows(mkIndex, mkIndex) const: subscript indices must be positive!" << endl;
        Basics::quit(1);
    }
    if (i1 > nrows || i2 > nrows) {
        *Basics::err << "SymmetricMatrix::columns(mkIndex, mkIndex) const: index exceeds matrix dimensions!" << endl;
        Basics::quit(1);
    }

    // find dimensions
    mkIndex nr = (i1<=i2) ? (i2-i1+1u) : (i1-i2+1u);

    // output storage address
    Real *sb = NULL;

    // copy elements
    submatrix_of_symmatrix(nrows, store, i1, i2, 1, nrows, sb);

    // return the result
    return Matrix(nr, nrows, sb);
}

// a submatrix formed by rows i1,...,i2 of *this
// this member function invokes submatrix_of_symmatrix()
Matrix SymmetricMatrix::subMatrix(mkIndex i1, mkIndex i2, mkIndex j1, mkIndex j2) const {
    if (i1 == 0 || i2 == 0 || j1 == 0 || j2 == 0) {
        *Basics::err << "SymmetricMatrix::subMatrix(mkIndex, mkIndex, mkIndex, mkIndex) const: subscript indices must be positive!" << endl;
        Basics::quit(1);
    }
    if (i1 > nrows || i2 > nrows || j1 > nrows || j2 > nrows) {
        *Basics::err << "SymmetricMatrix::subMatrix(mkIndex, mkIndex, mkIndex, mkIndex) const: index exceeds matrix dimensions!" << endl;
        Basics::quit(1);
    }

    // find dimensions
    mkIndex nr = (i1<=i2) ? (i2-i1+1u) : (i1-i2+1u);
    mkIndex nc = (j1<=j2) ? (j2-j1+1u) : (j1-j2+1u);

    // output storage address
    Real *sb = NULL;

    // copy elements
    submatrix_of_symmatrix(nrows, store, i1, i2, j1, j2, sb);

    // return the result
    return Matrix(nr, nc, sb);
}

// a submatrix as a symmetric matrix formed by elements in rows k1,...,k2 and columns k1,...,k2 of *this
// this member function invokes submatrix_of_symmatrix()
SymmetricMatrix SymmetricMatrix::subMatrix(mkIndex k1, mkIndex k2) const {
    // dimension of the submatrix
    mkIndex nrc = (k1<=k2) ? (k2-k1+1u) : (k1-k2+1u);

    Real *sb = NULL;
    submatrix_of_symmatrix(nrows, store, k1, k2, sb);
    return SymmetricMatrix(nrc, sb);
}


// return a vector formed by the elements in the k-th diagonal of *this in order
// this member function invokes symmatrix_diag()
Vector SymmetricMatrix::diag(mkSignedIndex k) const {
    // pointer for output
    Real *sv = NULL;

    // extract the elements in the k-th diagonal of the *this
    mkIndex len = symmatrix_diag(nrows, store, sv, k);

    // return the result
    return Vector(len, sv);
}


// Frobenius norm of *this
Real SymmetricMatrix::normFrobenius() const {
    return sqrt(sumSquare());
}

// sum of squared elements of *this
Real SymmetricMatrix::sumSquare() const {
    Real sum = 0.0;
    Real *sa = store;
    mkIndex nn = nrows;
    while (nn--) {
        // add up the diagonal element A(i,i), where i=nrows-nn
        sum += Basics::square(*sa++);
        mkIndex len = nn;
        while (len--) {
            // add up the off-diagonal element A(j,i) and also A(i,j), where j=i+(nn-len)
            sum += 2*Basics::square(*sa++);
        }
    }
    return sum;
}

// 1-norm of *this
Real SymmetricMatrix::norm1() const {
    // allocate work space
    Real *work = new Real[nrows];

    // zero initialization
    Real *w = work;
    mkIndex ii = nrows;
    while (ii--)
        *w++ = 0.0;

    // compute work[j] as the 1-norm of (j+1)st column
    Real *sa = store;
    for (mkIndex j=0; j<nrows; j++) {
        w = work+j;
        ii = nrows-j;
        while (ii--)
            *w++ += Basics::abs(*sa++);
    }

    // the 1-norm of *this is the maximum work[i] for i=0,...,nrows-1
    Real nrm = 0.0;
    for (mkIndex i=0; i<nrows; i++)
        nrm = (nrm>=work[i] ? nrm : work[i]);

    // free work space
    delete [] work;

    // return the result
    return nrm;
}

// infinity-norm of *this
Real SymmetricMatrix::normInfinity() const {
    return norm1();
}


// access A(i,j)
// return a reference to A(i,j), i.e. can be used to write A(i,j), where A is *this
Real& SymmetricMatrix::operator()(mkIndex i, mkIndex j) {
    if (i==0 || j==0) {
        *Basics::err << "SymmetricMatrix::operator()(mkIndex, mkIndex): subscript indices must be positive!" << endl;
        Basics::quit(1);
    }
    if (i>nrows || j>nrows) {
        *Basics::err << "SymmetricMatrix::operator()(mkIndex, mkIndex): index exceeds matrix dimensions!" << endl;
        Basics::quit(1);
    }

    // return a reference to A(i,j) if i>=j or A(j,i) if i<j, where A is *this
    if (i >= j)
        return store[(j-1u)*(2u*nrows-j+2u)/2u+(i-j)];
    return store[(i-1u)*(2*nrows-i+2u)/2u+(j-i)];  // write A(m,n)
}

// read A(i,j) via a const reference or pointer to A, where A is *this
Real SymmetricMatrix::operator()(mkIndex i, mkIndex j) const {
    if (i==0 || j==0) {
        *Basics::err << "SymmetricMatrix::operator()(mkIndex, mkIndex) const: subscript indices must be positive!" << endl;
        Basics::quit(1);
    }
    if (i>nrows || j>nrows) {
        *Basics::err << "SymmetricMatrix::operator()(mkIndex, mkIndex) const: index exceeds matrix dimensions!" << endl;
        Basics::quit(1);
    }

    // read A(i,j), via a const pointer or reference to A
    if (i >= j)
        return store[(j-1u)*(2u*nrows-j+2u)/2u+(i-j)];
    return store[(i-1u)*(2u*nrows-i+2u)/2u+(j-i)];  // read A(m,n), via a const reference or pointer
}

// operator =
// return a const reference to A after setting A(i,j)=f for i,j=1,...,nrows, where A is *this
const SymmetricMatrix& SymmetricMatrix::operator=(Real f) {
    // size of the storage to store the symmetric matrix *this
    mkIndex sz = nrows*(nrows+1u)/2u;

    // set A(i,j)=f for 1<=j<=i<=nrows
    Real *s = store;
    while (sz--)
        *(s++) = f;

    // return the result
    return *this;
}

// return a const reference to A after setting A=B, where A is *this
// this operator invokes memory_xcopy()
const SymmetricMatrix& SymmetricMatrix::operator=(const SymmetricMatrix &B) {
    // memory size
    mkIndex sz = B.Nrows()*(B.Nrows()+1u)/2u;

    if (nrows != B.Nrows()) {
        // free memory and allocate memory precisely required
        delete [] store;
        if (sz)
            store = new Real[sz];
        else
            store = NULL;
        nrows = B.Nrows();
    }

    // copy elements
    memory_xcopy(sz, B.Store(), store);

    // return the result
    return *this;
}

// operator -
// return -A, where A is *this
// this operator invokes vector_scalar_operation()
SymmetricMatrix SymmetricMatrix::operator-() const {
    // pointer for output
    Real *sb = NULL;

    // elementwise multiply -1.0
    vector_scalar_operation(nrows*(nrows+1u)/2u, Store(), -1.0, sb, 2);

    // return the result
    return SymmetricMatrix(nrows, sb);
}

// operators +,-,*,/,+=,-=,*-,/* a real number
// return B with B(i,j)=A(i,j)+f for i,j=1,...,nrows, where A is *this
// this operator invokes vector_scalar_operation()
SymmetricMatrix SymmetricMatrix::operator+(Real f) const {
    // pointer for output
    Real *sb = NULL;

    // compute B = A + f (elementwise), where A is *this and B is stored in sb[]
    vector_scalar_operation(nrows*(nrows+1u)/2u, Store(), f, sb, 0);

    // return the result
    return SymmetricMatrix(nrows, sb);
}

// return B with B(i,j)=A(i,j)-f for i,j=1,...,nrows, where A is *this
// this operator invokes vector_scalar_operation()
SymmetricMatrix SymmetricMatrix::operator-(Real f) const {
    // pointer for output
    Real *sb = NULL;

    // compute B = A - f (elementwise), where A is *this and B is stored in sb[]
    vector_scalar_operation(nrows*(nrows+1u)/2u, Store(), f, sb, 1);

    // return the result
    return SymmetricMatrix(nrows, sb);
}

// return B with B(i,j)=A(i,j)*f for i,j=1,...,nrows, where A is *this
// this operator invokes vector_scalar_operation()
SymmetricMatrix SymmetricMatrix::operator*(Real f) const {
    // pointer for output
    Real *sb = NULL;

    // compute B = A * f (elementwise), where A is *this and B is stored in sb[]
    vector_scalar_operation(nrows*(nrows+1u)/2u, Store(), f, sb, 2);

    // return the result
    return SymmetricMatrix(nrows, sb);
}

// return B with B(i,j)=A(i,j)/f for i,j=1,...,nrows, where A is *this
// this operator invokes vector_scalar_operation()
SymmetricMatrix SymmetricMatrix::operator/(Real f) const {
    // exception handling
    if (f == 0.0) {
        *Basics::err << "SymmetricMatrix::operator/(Real) const: division-by-zero exception!" << endl;
        Basics::quit(1);
    }

    // pointer for output
    Real *sb = NULL;

    // compute B = A / f (elementwise), where A is *this and B is stored in sb[]
    vector_scalar_operation(nrows*(nrows+1u)/2u, Store(), f, sb, 3);

    // return the result
    return SymmetricMatrix(nrows, sb);
}

// return a const reference to A after computing A(i,j)+=f for i,j=1,...,nrows, where A is *this
// this operator invokes vector_scalar_operation()
const SymmetricMatrix& SymmetricMatrix::operator+=(Real f) {
    // compute A += f (elementwise), where A is *this
    vector_scalar_operation(nrows*(nrows+1u)/2u, Store(), f, store, 0);

    // return the result
    return *this;
}

// return a const reference to A after computing A(i,j)-=f for i,j=1,...,nrows, where A is *this
// this operator (eventually) invokes vector_scalar_operation()
const SymmetricMatrix& SymmetricMatrix::operator-=(Real f) {
    return this->operator+=(-f);
}

// return a const reference to A after computing A(i,j)*=f for i,j=1,...,nrows, where A is *this
// this operator invokes vector_scalar_operation()
const SymmetricMatrix& SymmetricMatrix::operator*=(Real f) {
    // compute A *= f (elementwise), where A is *this
    vector_scalar_operation(nrows*(nrows+1u)/2u, Store(), f, store, 2);

    // return the result
    return *this;
}

// return a const reference to A after computing A(i,j)/=f for i,j=1,...,nrows, where A is *this
// this operator invokes vector_scalar_operation()
const SymmetricMatrix& SymmetricMatrix::operator/=(Real f) {
    // exception handling
    if (f == 0.0) {
        *Basics::err << "SymmetricMatrix::operator/=(Real): division-by-zero exception!" << endl;
        Basics::quit(1);
    }

    // compute A /= f (elementwise), where A is *this
    vector_scalar_operation(nrows*(nrows+1u)/2u, Store(), f, store, 3);

    // return the result
    return *this;
}

// scalar - symmetric operations
// return f*A
// this operator (eventually) invokes vector_scalar_operation()
SymmetricMatrix operator*(Real f, const SymmetricMatrix &A) {
    return A*f;
}

// return B with B(i,j)=f+A(i,j) for i,j=1,...,A.nrows
// this operator (eventually) invokes vector_scalar_operation()
SymmetricMatrix operator+(Real f, const SymmetricMatrix &A) {
    return A+f;
}

// return B with B(i,j)=f-A(i,j) for i,j=1,...,A.nrows
// this operator (eventually) invokes vector_scalar_operation()
SymmetricMatrix operator-(Real f, const SymmetricMatrix &A) {
    return (-A)+f;
}


// symmetric matrix - symmetric matrix addition/subtraction
// return A+B, where A is *this
// this operator invokes elementwise_addition()
SymmetricMatrix SymmetricMatrix::operator+(const SymmetricMatrix &B) const {
    // exception handling
    if (nrows != B.Nrows()) {
        *Basics::err << "SymmetricMatrix::operator+(const SymmetricMatrix&) const: matrix dimensions must agree!" << endl;
        Basics::quit(1);
    }

    // pointer for output
    Real *sc = NULL;

    // compute C = A + B, where A is *this and C stored in sc[]
    elementwise_addition(nrows*(nrows+1u)/2u, Store(), 1.0, B.Store(), sc);

    // return the result
    return SymmetricMatrix(nrows, sc);
}

// return A-B, where A is *this
// this operator invokes elementwise_addition()
SymmetricMatrix SymmetricMatrix::operator-(const SymmetricMatrix &B) const {
    // exception handling
    if (nrows != B.Nrows()) {
        *Basics::err << "SymmetricMatrix::operator-(const SymmetricMatrix&) const: matrix dimensions must agree!" << endl;
        Basics::quit(1);
    }

    // pointer for output
    Real *sc = NULL;

    // compute C = A - B, where A is *this and C is stored in sc[]
    elementwise_addition(nrows*(nrows+1u)/2u, Store(), -1.0, B.Store(), sc);

    // return the result
    return SymmetricMatrix(nrows, sc);
}

// return a const reference to A after computing A = A+B, where A is *this
// this operator invokes elementwise_addition()
const SymmetricMatrix& SymmetricMatrix::operator+=(const SymmetricMatrix &B) {
    if (nrows != B.Nrows()) {
        *Basics::err << "SymmetricMatrix::operator+=(const SymmetricMatrix&): matrix dimensions must agree!" << endl;
        Basics::quit(1);
    }

    // compute A += B, where A is *this
    elementwise_addition(nrows*(nrows+1u)/2u, Store(), 1.0, B.Store(), store);

    // return the result
    return *this;
}

// return a const reference to A after computing A = A-B, where A is *this
// this operator invokes elementwise_addition()
const SymmetricMatrix& SymmetricMatrix::operator-=(const SymmetricMatrix &B) {
    // exception handling
    if (nrows != B.Nrows()) {
        *Basics::err << "SymmetricMatrix::operator-=(const SymmetricMatrix&): matrix dimensions must agree!" << endl;
        Basics::quit(1);
    }

    // compute A -= B, where A is *this
    elementwise_addition(nrows*(nrows+1u)/2u, Store(), -1.0, B.Store(), store);

    // return the result
    return *this;
}

// symmetric matrix - matrix addition/subtraction
// return A+B, where A is *this
// this operator (eventually) invokes matrix_symmatrix_addition()
Matrix SymmetricMatrix::operator+(const Matrix &B) const {
    if (nrows!=B.Nrows() || nrows!=B.Ncols()) {
        *Basics::err << "SymmetricMatrix::operator+(const Matrix&) const: matrix dimensions must agree!" << endl;
        Basics::quit(1);
    }

    return B+(*this);
}

// return A-B, where A is *this
// this operator (eventually) invokes matrix_symmatrix_addition()
Matrix SymmetricMatrix::operator-(const Matrix &B) const {
    if (nrows!=B.Nrows() || nrows!=B.Ncols()) {
        *Basics::err << "SymmetricMatrix::operator-(const Matrix&) const: matrix dimensions must agree!" << endl;
        Basics::quit(1);
    }

    return -B+(*this);
}


// symmetric matrix - vector multiplication
// return A*v, where A is *this
// this operator invokes symmatrix_vector_multiplication()
Vector SymmetricMatrix::operator*(const Vector &v) const {
    if (nrows != v.Length()) {
        *Basics::err << "SymmetricMatrix::operator*(const Vector&) const: inner matrix-vector dimensions must agree!" << endl;
        Basics::quit(1);
    }
#ifdef USE_MWBLAS
    // MWBLAS routines sspmv and dspmv require non-const storage for reading the input matrix A and vector v,
    // so it cannot be incorporated into matrix_vector_multiplication and therefore is put here
    Real *sw = new Real[nrows];
    char uplo = 'L';
    Real alp = 1.0, bet = 0.0;
    mwSignedIndex inc = 1;
    mwSignedIndex nrows2 = (mwSignedIndex)nrows;
    #ifdef USE_SINGLE
        sspmv(&uplo, &nrows2, &alp, store, v.store, &inc, &bet, sw, &inc);
    #else
        dspmv(&uplo, &nrows2, &alp, store, v.store, &inc, &bet, sw, &inc);
    #endif
#else
    // pointer for output
    Real *sw = NULL;

    // compute w = A*v, where A is *this and w is stored in sw[]
    symmatrix_vector_multiplication(nrows, store, v.Store(), sw);
#endif

    // return the result
    return Vector(nrows, sw);
}


// symmetric matrix - symmetric matrix multiplication
// return A*B, where A is *this
// this operator invokes symmatrix_symmatrix_multiplication()
Matrix SymmetricMatrix::operator*(const SymmetricMatrix &B) const {
    if (nrows != B.Nrows()) {
        *Basics::err << "SymmetricMatrix::operator*(const SymmetricMatrix&) const: inner matrix dimensions must agree!" << endl;
        Basics::quit(1);
    }

    // pointer for output
    Real *sc = NULL;

    // compute C = A*B, where A is *this and C is stored in sc[]
    symmatrix_symmatrix_multiplication(nrows, Store(), B.Store(), sc);

    // return the result
    return Matrix(nrows, nrows, sc);
}

// symmetric matrix - matrix multiplication
// return A*B, where A is *this
// this operator invokes symmatrix_matrix_multiplication()
Matrix SymmetricMatrix::operator*(const Matrix &B) const {
    // // exception handling
    if (nrows != B.Nrows()) {
        *Basics::err << "SymmetricMatrix::operator*(const Matrix&) const: inner matrix dimensions must agree!" << endl;
        Basics::quit(1);
    }

    // pointer for output
    Real *sc = NULL;

    // compute C = A*B, where A is *this and C is stored in sc[]
    symmatrix_matrix_multiplication(nrows, B.Ncols(), Store(), B.Store(), sc);

    // return the result
    return Matrix(nrows, B.Ncols(), sc);
}


// symmetrically permute rows and columns of A, such that
//    B(perm[i-1],perm[j-1])     = A(i,j) for i,j=1,...,nrows (if indexFrom1==true), or
//    B(perm[i-1]+1,perm[j-1]+1) = A(i,j) for i,j=1,...,nrows (if indexFrom1==false),
// where A is *this and B is the return matrix
// the indices perm[0,...,nrows-1] range over 0,1,...,nrows-1 (if indexFrom1==false) or over 1,2,...,nrows (if indexFrom1==true)
// this member function invokes permute_symmatrix_rows_and_columns()
// this member function resembles B(perm,perm) = A in OCTAVE/MATLAB
SymmetricMatrix SymmetricMatrix::permuteRowsAndColumns(const mkIndex *perm, bool indexFrom1) const {
    // pointer for output
    Real *sb = NULL;

    // permute
    permute_symmatrix_rows_and_columns(nrows, perm, store, sb, indexFrom1);

    // return the result
    return SymmetricMatrix(nrows, sb);
}


#ifdef USE_NAMESPACE
}  // end of namespace NEWMAT
#endif
