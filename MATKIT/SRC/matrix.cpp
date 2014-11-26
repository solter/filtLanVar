//============================================================================
// the MATKIT class Matrix
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


Real Matrix::memoryExpansionFactor = 1.2;


// display matrix A with operator <<
// in other words, print out matrix A to ostream s
// this routine invokes print_matrix()
std::ostream& operator<<(std::ostream &s, const Matrix &A) {
    print_matrix(s, A.Nrows(), A.Ncols(), A.Store());
    return s;
}

// return an nr-by-nc matrix of zeros
// this static function resembles zeros(nr,nc) in OCTAVE/MATLAB
Matrix Matrix::zeros(mkIndex nr, mkIndex nc) {
    // allocate memory
    mkIndex sz = nr*nc;
    Real *sa = new Real[sz];

    // set it zero
    Real *a0 = sa;
    // memset(sa, 0, sz);  // works for IEEE 754, not guaranteed working for other types of floating point numbers
    while (sz--)
        *a0++ = 0.0;

    // return the result
    return Matrix(nr, nc, sa);
}

// return an nr-by-nc matrix of ones
// this static function resembles ones(nr,nc) in OCTAVE/MATLAB
Matrix Matrix::ones(mkIndex nr, mkIndex nc) {
    // allocate memory
    mkIndex sz = nr*nc;
    Real *sa = new Real[sz];

    // set ones
    Real *a0 = sa;
    while (sz--)
        *a0++ = 1.0;

    // return the result
    return Matrix(nr, nc, sa);
}

// return an nr-by-nc matrix with ones on the main diagonal and zeros elsewhere
// this static function resembles eye(nr,nc) in OCTAVE/MATLAB
Matrix Matrix::eye(mkIndex nr, mkIndex nc) {
    // allocate memory
    mkIndex sz = nr*nc;
    Real *sa = new Real[sz];

    // zero initialization
    Real *a0 = sa;
    while (sz--)
        *a0++ = 0.0;

    // set the eye
    a0 = sa;
    mkIndex step = nr+1u;
    mkIndex mrc = (nr<=nc) ? nr : nc;
    while (mrc--) {
        *a0 = 1.0;
        a0 += step;
    }

    // return the result
    return Matrix(nr, nc, sa);
}

// return an nr-by-nc matrix containing pseudo-random values drawn from a uniform distribution in the interval [a,b]
// this static function resembles rand(nr,nc) in OCTAVE/MATLAB
Matrix Matrix::random(mkIndex nr, mkIndex nc, Real a, Real b) {
    // allocate memory
    mkIndex sz = nr*nc;
    Real *sa = new Real[sz];

    // set random elements
    Real *a0 = sa;
    while (sz--)
        *a0++ = a+(b-a)*((Real)rand()/(Real)RAND_MAX);

    // return the result
    return Matrix(nr, nc, sa);
}


// horizontal concatenation of two matrices A and B
// A and B must have the same number of rows
// this static function invokes horizontal_concatenate()
// this static function resembles [A,B] or equivalently horzcat(A,B) in OCTAVE/MATLAB
Matrix Matrix::horizontalConcatenate(const Matrix &A, const Matrix &B) {
    // exception handling
    mkIndex nr = A.Nrows();
    if (B.Nrows() != nr) {
        // the two matrices should have the same number of rows
        *Basics::err << "Matrix horizontalConcatenate(const Matrix A&, const Matrix &B): A and B should have the same number of rows!" << endl;
        Basics::quit(1);
    }

    // pointer for the output matrix
    Real *sc = NULL;

    // copy elements
    horizontal_concatenate(nr, A.Ncols(), A.Store(), B.Ncols(), B.Store(), sc);

    // return the result
    return Matrix(nr, A.Ncols()+B.Ncols(), sc);
}

// vertical concatenation of two matrices A and B
// A and B must have the same number of columns
// this static function invokes vertical_concatenate()
// this static function resembles [A;B] or equivalently vertcat(A,B) in OCTAVE/MATLAB
Matrix Matrix::verticalConcatenate(const Matrix &A, const Matrix &B) {
    // exception handling
    mkIndex nc = A.Ncols();
    if (B.Ncols() != nc) {
        // the two matrices should have the same number of rows
        *Basics::err << "Matrix verticalConcatenate(const Matrix A&, const Matrix &B): A and B should have the same number of columns!" << endl;
        Basics::quit(1);
    }

    // pointer for the output matrix
    Real *sc = NULL;

    // copy elements
    vertical_concatenate(A.Nrows(), nc, A.Store(), B.Nrows(), B.Store(), sc);

    // return the result
    return Matrix(A.Nrows()+B.Nrows(), nc, sc);
}

// sum of two matrices, padding zeros if they have different dimensions
// this static function invokes matrix_xsum()
Matrix Matrix::xsum(const Matrix &A, const Matrix &B) {
    // find dimensions
    mkIndex nr = Basics::max(A.Nrows(), B.Nrows());
    mkIndex nc = Basics::max(A.Ncols(), B.Ncols());

    // pointer for output
    Real *sc = NULL;

    // compute xsum
    matrix_xsum(A.Nrows(), A.Ncols(), A.Store(), B.Nrows(), B.Ncols(), B.Store(), sc);

    // return the result
    return Matrix(nr, nc, sc);
}


// constructors
// a constructor for an empty matrix
Matrix::Matrix() {
    // a zero matrix
    nrows = 0;
    ncols = 0;
    store = NULL;
    maxncols = 0;
}

// a constructor for an nr-by-nc matrix of zeros
// if maxnc<=nc, then memory of size nr*nc will be allocated
// if maxnc>nc,  then memory of size nr*maxnc will be allocated
Matrix::Matrix(mkIndex nr, mkIndex nc, mkIndex maxnc) {
    // set dimensions
    nrows = nr;
    ncols = nc;
    maxncols = (maxnc>=nc) ? maxnc : nc;

    // allocate memory
    if (nr > 0 && maxncols > 0)
        store = new Real[nr*maxncols];
    else
        store = NULL;

    // zero initialization
    if (nr > 0 && nc > 0) {
        mkIndex sz = nr*nc;
        Real *sa = store;
        while (sz--)
            *sa++ = 0.0;
    }
}

// a constructor for an nr-by-nc matrix of with elements stored in store0[]
// if maxnc<=nc, then it is assumed that memory of size at least nr*nc    has been allocated, with the address pointed to by store0
// if maxnc>nc,  then it is assumed that memory of size at least nr*maxnc has been allocated, with the address pointed to by store0
// note that it performs a shallow copy of store0 and the memory will be freed by the destructor (i.e. `delete [] store0')
Matrix::Matrix(mkIndex nr, mkIndex nc, Real *store0, mkIndex maxnc) {
    // set dimensions
    nrows = nr;
    ncols = nc;

    // a shallow copy
    store = store0;

    // set maximum number of columns allowed in the memory pointed to by store0
    maxncols = (maxnc>=nc) ? maxnc : nc;
}

// a copy constructor
// this constructor invokes memory_xcopy()
Matrix::Matrix(const Matrix &B) {
    // set dimensions
    nrows = B.Nrows();
    ncols = B.Ncols();
    maxncols = ncols;

    // allocate memory
    mkIndex sz = nrows*ncols;
    store = new Real[sz];

    // copy the elements
    memory_xcopy(sz, B.Store(), store);
}

// a destructor
Matrix::~Matrix() {
    // free memory
    delete [] store;
}


// convert *this to a sparse matrix
// this member function invokes full_to_sparse()
// this member function resembles sparse(A) in OCTAVE/MATLAB, where A is *this
SparseMatrix Matrix::toSparse() const {
    // pointers for output
    Real *sa = NULL;
    mkIndex *rowIdx = NULL;
    mkIndex *colPtr = NULL;

    // convert *this to a sparse matrix stored in sa[], rowIdx[], colPtr[]
    full_to_sparse(nrows, ncols, Store(), sa, rowIdx, colPtr);

    // return the result
    return SparseMatrix(nrows, ncols, sa, rowIdx, colPtr);
}

// resize *this to be of size nr-by-nc
// if lazyMalloc == false, then memory will be allocated whenever the size of allocated memory is different from required
// if lazyMalloc == true, then memory will be allocated only when the allocated memory is insufficient
void Matrix::resize(mkIndex nr, mkIndex nc, bool lazyMalloc) {
    // memory size already allocated
    mkIndex sz0 = nrows*maxncols;
    // memory size required
    mkIndex sz = nr*nc;

    if (lazyMalloc) {
        // allocate memory only if it is insufficient
        if (sz > sz0) {
            // free memory, and allocate memory precisely required
            delete [] store;
            if (sz > 0)
                store = new Real[sz];
            else
                store = NULL;
            maxncols = nc;
        }
        else {
            maxncols = sz0 / nr;
        }
    }
    else {
        // use memory of precisely required size
        if (sz != sz0) {
            // free memory, and allocate memory precisely required
            delete [] store;
            if (sz > 0)
                store = new Real[sz];
            else
                store = NULL;
        }
        // set maximum number of columns allowed in the allocated memory
        maxncols = nc;
    }

    // set dimensions
    nrows = nr;
    ncols = nc;
}

// set the matrix
// if maxnc<=nc, then it is assumed that memory of size at least nr*nc    has been allocated, with the address pointed to by store0
// if maxnc>nc,  then it is assumed that memory of size at least nr*maxnc has been allocated, with the address pointed to by store0
// it performs a shallow copy of store0[] and the memory will be freed by the destructor (i.e. `delete [] store0')
void Matrix::set(mkIndex nr, mkIndex nc, Real *store0, mkIndex maxnc) {
    // set dimensions
    nrows = nr;
    ncols = nc;
    maxncols = (maxnc>=nc) ? maxnc : nc;

    if (store != store0) {
        // free memory, and do a shallow copy
        delete [] store;
        store = store0;
    }
}

// append a column formed by sv[0,...,nrows-1] to *this
// allocated memory, if insufficient, will be expanded by a factor of Matrix::memoryExpansionFactor
// this member function invokes memory_xcopy()
// see also Matrix::setMemoryExpansionFactor()
void Matrix::appendColumn(const Real *sv) {
    if (ncols == maxncols) {
        // insufficient memory, expand
        maxncols *= memoryExpansionFactor;
        if (maxncols <= ncols)
            maxncols ++;
        // allocate more memory
        Real *store2 = new Real[nrows*maxncols];
        // copy old elements
        memory_xcopy(nrows*ncols, Store(), store2);
        // free old memory and set the pointer
        delete [] store;
        store = store2;
    }

    // append the column
    memory_xcopy(nrows, sv, store+nrows*ncols);

    // now we have one more column
    ncols ++;
}

// append a column formed by v to *this
// allocated memory, if insufficient, will be expanded by a factor of Matrix::memoryExpansionFactor
// see also Matrix::setMemoryExpansionFactor()
void Matrix::appendColumn(const Vector &v) {
    // exception handling
    if (nrows != v.Length()) {
        *Basics::err << "Matrix::appendColumn(const Vector&): dimensions are not consistent!" << endl;
        Basics::quit(1);
    }

    // append vector v to be the last column of *this
    appendColumn(v.Store());
}

// static function setMemoryExpansionFactor
// the default memoryExpansionFactor is 1.2, shown in matrix.cpp
// the memoryExpansionFactor is used by Matrix::appendColumn()
void Matrix::setMemoryExpansionFactor(Real newExpansionFactor) {
    if (newExpansionFactor <= 1.0 || newExpansionFactor > 2.0) {
        *Basics::err << "Matrix::setMemoryExpansionFactor(Real): expansion factor should be >1.0 and <= 2.0!" << endl;
        Basics::quit(1);
    }

    // set the memory expansion factor
    memoryExpansionFactor = newExpansionFactor;
}

// a vector formed by column j of *this
// this member function invokes memory_xcopy()
// this member function resembles A(:,j) in OCTAVE/MATLAB, where A is *this
Vector Matrix::column(mkIndex j) const {
    // exception handling
    if (j == 0) {
        *Basics::err << "Matrix::column(mkIndex) const: subscript indices must be positive!" << endl;
        Basics::quit(1);
    }
    if (j > ncols) {
        *Basics::err << "Matrix::column(mkIndex) const: index exceeds matrix dimensions!" << endl;
        Basics::quit(1);
    }

    // allocate memory
    Real *sc = new Real[nrows];

    // copy elements
    memory_xcopy(nrows, Store()+(j-1u)*nrows, sc);

    // return the result
    return Vector(nrows, sc);
}

// a vector formed by row i of *this
// this member function invokes memory_xcopy()
// this member function resembles A(i,:) in OCTAVE/MATLAB, where A is *this
Vector Matrix::row(mkIndex i) const {
    // exception handling
    if (i == 0) {
        *Basics::err << "Matrix::row(mkIndex) const: subscript indices must be positive!" << endl;
        Basics::quit(1);
    }
    if (i > nrows) {
        *Basics::err << "Matrix::row(mkIndex) const: index exceeds matrix dimensions!" << endl;
        Basics::quit(1);
    }

    // allocate memory
    Real *sr = new Real[ncols];

    // copy elements
    memory_xcopy(ncols, Store()+(i-1u), sr, (mkSignedIndex)nrows, 1);

    // return the result
    return Vector(ncols, sr);
}

// a submatrix formed by columns j1,...,j2 of *this
// this member function invokes submatrix_of_matrix()
// this member function resembles A(:,j1:j2) in OCTAVE/MATLAB, where A is *this
Matrix Matrix::columns(mkIndex j1, mkIndex j2) const {
    // exception handling
    if (j1 == 0 || j2 == 0) {
        *Basics::err << "Matrix::columns(mkIndex, mkIndex) const: subscript indices must be positive!" << endl;
        Basics::quit(1);
    }
    if (j1 > ncols || j2 > ncols) {
        *Basics::err << "Matrix::columns(mkIndex, mkIndex) const: index exceeds matrix dimensions!" << endl;
        Basics::quit(1);
    }

    // find the number of columns of the submatrix
    mkIndex ncb = (j1<=j2) ? (j2-j1+1u) : (j1-j2+1u);

    // compute the submatrix
    Real *sb = NULL;
    submatrix_of_matrix(nrows, ncols, Store(), 1, nrows, j1, j2, sb);

    // return the result
    return Matrix(nrows, ncb, sb);
}

// a submatrix formed by rows i1,...,i2 of *this
// this member function invokes submatrix_of_matrix()
// this member function resembles A(i1:i2,:) in OCTAVE/MATLAB, where A is *this
Matrix Matrix::rows(mkIndex i1, mkIndex i2) const {
    // exception handling
    if (i1 == 0 || i2 == 0) {
        *Basics::err << "Matrix::rows(mkIndex, mkIndex) const: subscript indices must be positive!" << endl;
        Basics::quit(1);
    }
    if (i1 > nrows || i2 > nrows) {
        *Basics::err << "Matrix::rows(mkIndex, mkIndex) const: index exceeds matrix dimensions!" << endl;
        Basics::quit(1);
    }

    // find the number of rows of the submatrix
    mkIndex nrb = (i1<=i2) ? (i2-i1+1u) : (i1-i2+1u);

    // compute the submatrix
    Real *sb = NULL;
    submatrix_of_matrix(nrows, ncols, Store(), i1, i2, 1, ncols, sb);

    // return the result
    return Matrix(nrb, ncols, sb);
}

// a submatrix formed by elements in rows i1,...,i2 and columns j1,...,j2 of *this
// this member function invokes submatrix_of_matrix()
// this member function resembles A(i1:i2,j1:j2) in OCTAVE/MATLAB, where A is *this
Matrix Matrix::subMatrix(mkIndex i1, mkIndex i2, mkIndex j1, mkIndex j2) const {
    // exception handling
    if (i1 == 0 || i2 == 0 || j1 == 0 || j2 == 0) {
        *Basics::err << "Matrix::subMatrix(mkIndex, mkIndex, mkIndex, mkIndex) const: subscript indices must be positive!" << endl;
        Basics::quit(1);
    }
    if (i1 > nrows || i2 > nrows || j1 > ncols || j2 > ncols) {
        *Basics::err << "Matrix::subMatrix(mkIndex, mkIndex, mkIndex, mkIndex) const: index exceeds matrix dimensions!" << endl;
        Basics::quit(1);
    }

    // find dimensions
    mkIndex nrb = (i1<=i2) ? (i2-i1+1u) : (i1-i2+1u);
    mkIndex ncb = (j1<=j2) ? (j2-j1+1u) : (j1-j2+1u);

    // compute the submatrix
    Real *sb = NULL;
    submatrix_of_matrix(nrows, ncols, Store(), i1, i2, j1, j2, sb);

    // return the result
    return Matrix(nrb, ncb, sb);
}

// a submatrix formed by columns c[0],c[1],...,c[len-1] of *this
// this member function invokes memory_xcopy()
// this member function resembles A(:,c), where A is *this
Matrix Matrix::columns(mkIndex len, const mkIndex c[], bool indexFrom1) const {
    // allocate memory
    Real *sb = new Real[nrows*len];

    const Real *sa = store;
    if (indexFrom1) {
        // pre-shift the pointer
        sa -= nrows;
    }

    // compute B = A(:,c[0,...,len-1])
    for (mkIndex jj=0; jj<len; jj++) {
        mkIndex j = c[jj];
        // exception handling
        if (indexFrom1 && j == 0) {
            *Basics::err << "Matrix::columns(mkIndex, const mkIndex*, bool indexFrom1) const: subscript indices must be positive if indexFrom1==true!" << endl;
            Basics::quit(1);
        }
        if (j >= ncols+indexFrom1) {
            *Basics::err << "Matrix::columns(mkIndex, const mkIndex*, bool) const: index exceeds matrix dimensions!" << endl;
            Basics::quit(1);
        }

        // compute B(:,jj+1) = A(:,j)
        const Real *a0 = sa + j*nrows;
        Real *b0 = sb + jj*nrows;
        memory_xcopy(nrows, a0, b0);
    }

    // return the result
    return Matrix(nrows, len, sb);
}

// a submatrix formed by rows r[0],r[1],...,r[len-1] of *this
// this member function resembles A(r,:), where A is *this
Matrix Matrix::rows(mkIndex len, const mkIndex r[], bool indexFrom1) const {
    // check the indices first
    mkIndex ii = len;
    for (mkIndex ii=0; ii<len; ii++) {
        if (indexFrom1 && r[ii] == 0) {
            *Basics::err << "Matrix::rows(mkIndex, const mkIndex*, bool indexFrom1) const: subscript indices must be positive if indexFrom1==true!" << endl;
            Basics::quit(1);
        }
        if (r[ii] >= nrows+indexFrom1) {
            *Basics::err << "Matrix::rows(mkIndex, const mkIndex*, bool) const: index exceeds matrix dimensions!" << endl;
            Basics::quit(1);
        }
    }

    // allocate memory
    Real *sb = new Real[len*ncols];

    const Real *sa = store;
    if (indexFrom1) {
        // pre-shift the pointer
        sa --;
    }

    // compute B = A(r[0,...,len-1],:)
    for (mkIndex j=0; j<ncols; j++) {
        // compute B(:,j+1) = A(r[0,...,len-1],j+1)
        const mkIndex *ri = r;
        const Real *a0 = sa + j*nrows;
        Real *b0 = sb + j*len;
        ii = len;
        while (ii--)
            *b0++ = a0[*ri++];
    }

    // return the result
    return Matrix(len, ncols, sb);
}

// a submatrix formed by elements in rows r[0],r[1],...,r[rlen-1] and columns c[0],c[1],...,c[clen-1] of *this
// this member function resembles A(r,c), where A is *this
Matrix Matrix::subMatrix(mkIndex rlen, const mkIndex r[], mkIndex clen, const mkIndex c[], bool indexFrom1) const {
    return rows(rlen,r,indexFrom1).columns(clen,c,indexFrom1);
}


// a vector formed by the elements in the k-th diagonal of *this in order
// this member function invokes matrix_diag()
// this member function resembles diag(A,k) in OCTAVE/MATLAB
Vector Matrix::diag(mkSignedIndex k) const {
    // pointer for output
    Real *sd = NULL;

    // extract the k-th diagonal
    mkIndex len = matrix_diag(nrows, ncols, store, sd, k);

    // return the result
    return Vector(len, sd);
}

// transpose of *this
// this member function invokes matrix_transpose()
// this member function resembles A' or equivalently transpose(A) in OCTAVE/MATLAB, where A is *this
Matrix Matrix::transpose() const {
    // pointer for the output
    Real *sb = NULL;

    // transpose *this
    matrix_transpose(nrows, ncols, Store(), sb);

    // return the result
    return Matrix(ncols, nrows, sb);
}


// Frobenius norm of *this
Real Matrix::normFrobenius() const {
#if defined(USE_MWBLAS)
    mwSignedIndex sz = (mwSignedIndex)(nrows*ncols);
    mwSignedIndex inc = 1;
    #ifdef USE_SINGLE
        return snrm2(&sz, store, &inc);
    #else
        return dnrm2(&sz, store, &inc);
    #endif
#elif defined(USE_BLAS)
    int sz = (int)(nrows*ncols);
    int inc = 1;
    #ifdef USE_SINGLE
        return snrm2_(&sz, store, &inc);
    #else
        return dnrm2_(&sz, store, &inc);
    #endif
#else
    return sqrt(sumSquare());
#endif
}

// sum of squared elements of *this
Real Matrix::sumSquare() const {
    // initialization
    Real sum = 0.0;

    // sum up squared elements
    Real *s = store;
    mkIndex sz = nrows*ncols;
    while (sz--)
        sum += Basics::square(*s++);

    // return the result
    return sum;
}

// 1-norm of *this
Real Matrix::norm1() const {
    // initialization
    Real nrm = 0.0;

    // compute the 1-norm of *this
    Real *sa = store;
    for (mkIndex j=1; j<=ncols; j++) {
        // compute the 1-norm of column j
        Real val = 0.0;
        mkIndex ii = nrows;
        while (ii--)
            val += Basics::abs(*sa++);
        // the matrix 1-norm is the maximum column 1-norm
        nrm = (nrm>=val) ? nrm : val;
    }

    // return the result
    return nrm;
}

// infinity-norm of *this
Real Matrix::normInfinity() const {
    // allocate work space of size nrows and initialize it zero
    Real *work = new Real[nrows];

    // zero initialization
    Real *w = work;
    mkIndex ii = nrows;
    while (ii--)
        *w++ = 0.0;

    // compute the infinity norm of each row
    Real *sa = store;
    for (mkIndex j=1; j<=ncols; j++) {
        w = work;
        ii = nrows;
        while (ii--)
            *w++ += Basics::abs(*sa++);
    }

    // the matrix infinity norm is the maximum row infinity norm
    Real nrm = 0.0;
    for (mkIndex i=0; i<nrows; i++)
        nrm = (nrm>=work[i]) ? nrm : work[i];

    // free the work space
    delete [] work;

    // return the result
    return nrm;
}


// access A(i,j)
// return a reference to A(i,j), i.e. can be used to write A(i,j), where A is *this
Real& Matrix::operator()(mkIndex i, mkIndex j) {
    // exception handling
    if (i==0 || j==0) {
        *Basics::err << "Matrix::operator()(mkIndex, mkIndex): subscript indices must be positive!" << endl;
        Basics::quit(1);
    }
    if (i>nrows || j>ncols) {
        *Basics::err << "Matrix::operator()(mkIndex, mkIndex): index exceeds matrix dimensions!" << endl;
        Basics::quit(1);
    }

    // return a reference to A(i,j), where A is *this
    return store[(i-1u)+(j-1u)*nrows];
}

// read A(i,j) via a const reference or pointer to A, where A is *this
Real Matrix::operator()(mkIndex i, mkIndex j) const {
    // exception handling
    if (i==0 || j==0) {
        *Basics::err << "Matrix::operator()(mkIndex, mkIndex) const: subscript indices must be positive!" << endl;
        Basics::quit(1);
    }
    if (i>nrows || j>ncols) {
        *Basics::err << "Matrix::operator()(mkIndex, mkIndex) const: index exceeds matrix dimensions!" << endl;
        Basics::quit(1);
    }

    // read A(i,j), via a const pointer or reference to A
    return store[(i-1u)+(j-1u)*nrows];
}


// operator =
// return a const reference to A after setting A(i,j)=f for i=1,...,nrows and j=1,...,ncols, where A is *this
const Matrix& Matrix::operator=(Real f) {
    // set all elements be f
    mkIndex sz = nrows*ncols;
    Real *sa = store;
    while (sz--)
        *sa++ = f;

    // return the result
    return *this;
}

// return a const reference to A after setting A=B, where A is *this
// this operator invokes memory_xcopy()
const Matrix& Matrix::operator=(const Matrix &B) {
    // size of the output matrix
    mkIndex sz = B.Nrows()*B.Ncols();

    if (nrows*ncols != sz) {
        // free memory and allocate memory precisely required
        delete [] store;
        if (sz)
            store = new Real[sz];
        else
            store = NULL;
    }

    // set dimensions and copy elements
    nrows = B.Nrows();
    ncols = B.Ncols();
    memory_xcopy(sz, B.Store(), store);

    // return the result
    return *this;
}

// operator -
// return -A, where A is *this
// this operator invokes vector_scalar_operation()
Matrix Matrix::operator-() const {
    // pointer for output
    Real *sb = NULL;

    // elementwise multiply -1.0
    vector_scalar_operation(nrows*ncols, Store(), -1.0, sb, 2);

    // return the result
    return Matrix(nrows, ncols, sb);
}

// operators +,-,*,/,+=,-=,*=,/= a real number
// return B with B(i,j)=A(i,j)+f for i=1,...,nrows and j=1,...,ncols, where A is *this
// this operator invokes vector_scalar_operation()
Matrix Matrix::operator+(Real f) const {
    // pointer for output
    Real *sb = NULL;

    // compute B = A + f (elementwise), where A is *this and B is stored in sb[]
    vector_scalar_operation(nrows*ncols, Store(), f, sb, 0);

    // return the result
    return Matrix(nrows, ncols, sb);
}

// return B with B(i,j)=A(i,j)-f for i=1,...,nrows and j=1,...,ncols, where A is *this
// this operator invokes vector_scalar_operation()
Matrix Matrix::operator-(Real f) const {
    // pointer for output
    Real *sb = NULL;

    // compute B = A - f (elementwise), where A is *this and B is stored in sb[]
    vector_scalar_operation(nrows*ncols, Store(), f, sb, 1);

    // return the result
    return Matrix(nrows, ncols, sb);
}

// return B with B(i,j)=A(i,j)*f for i=1,...,nrows and j=1,...,ncols, where A is *this
// this operator invokes vector_scalar_operation()
Matrix Matrix::operator*(Real f) const {
    // pointer for output
    Real *sb = NULL;

    // compute B = A * f (elementwise), where A is *this and B is stored in sb[]
    vector_scalar_operation(nrows*ncols, Store(), f, sb, 2);

    // return the result
    return Matrix(nrows, ncols, sb);
}

// return B with B(i,j)=A(i,j)/f for i=1,...,nrows and j=1,...,ncols, where A is *this
// this operator invokes vector_scalar_operation()
Matrix Matrix::operator/(Real f) const {
    // exception handling
    if (f == 0.0) {
        *Basics::err << "Matrix::operator/(Real) const: division-by-zero exception!" << endl;
        Basics::quit(1);
    }

    // pointer for output
    Real *sb = NULL;

    // compute B = A / f (elementwise), where A is *this and B is stored in sb[]
    vector_scalar_operation(nrows*ncols, Store(), f, sb, 3);

    // return the result
    return Matrix(nrows, ncols, sb);
}

// return a const reference to A after computing A(i,j)+=f for i=1,...,nrows and j=1,...,ncols, where A is *this
// this operator invokes vector_scalar_operation()
const Matrix& Matrix::operator+=(Real f) {
    // compute A += f (elementwise), where A is *this
    vector_scalar_operation(nrows*ncols, Store(), f, store, 0);

    // return the result
    return *this;
}

// return a const reference to A after computing A(i,j)-=f for i=1,...,nrows and j=1,...,ncols, where A is *this
// this operator (eventually) invokes vector_scalar_operation()
const Matrix& Matrix::operator-=(Real f) {
    return this->operator+=(-f);
}

// return a const reference to A after computing A(i,j)*=f for i=1,...,nrows and j=1,...,ncols, where A is *this
// this operator invokes vector_scalar_operation()
const Matrix& Matrix::operator*=(Real f) {
    // compute A *= f (elementwise), where A is *this
    vector_scalar_operation(nrows*ncols, Store(), f, store, 2);

    // return the result
    return *this;
}

// return a const reference to A after computing A(i,j)/=f for i=1,...,nrows and j=1,...,ncols, where A is *this
// this operator invokes vector_scalar_operation()
const Matrix& Matrix::operator/=(Real f) {
    // exception handling
    if (f == 0.0) {
        *Basics::err << "Matrix::operator/=(Real): division-by-zero exception!" << endl;
        Basics::quit(1);
    }

    // compute A /= f (elementwise), where A is *this
    vector_scalar_operation(nrows*ncols, Store(), f, store, 3);

    // return the result
    return *this;
}

// scalar - matrix operations
// return f*A
// this operator (eventually) invokes vector_scalar_operation()
Matrix operator*(Real f, const Matrix &A) {
    return A*f;
}

// return B with B(i,j)=f+A(i,j) for i=1,...,A.nrows and j=1,...,A.ncols
// this operator (eventually) invokes vector_scalar_operation()
Matrix operator+(Real f, const Matrix &A) {
    return A+f;
}

// return B with B(i,j)=f-A(i,j) for i=1,...,A.nrows and j=1,...,A.ncols
// this operator (eventually) invokes vector_scalar_operation()
Matrix operator-(Real f, const Matrix &A) {
    return (-A)+f;
}


// matrix - matrix addition/subtraction
// return A+B, where A is *this
// this operator invokes elementwise_addition()
Matrix Matrix::operator+(const Matrix &B) const {
    // exception handling
    if (nrows!=B.Nrows() || ncols!=B.Ncols()) {
        *Basics::err << "Matrix::operator+(const Matrix&) const: matrix dimensions must agree!" << endl;
        Basics::quit(1);
    }

    // pointer for output
    Real *sc = NULL;

    // compute C = A + B, where A is *this and C is stored in sc[]
    elementwise_addition(nrows*ncols, Store(), 1.0, B.Store(), sc);

    // return the result
    return Matrix(nrows, ncols, sc);
}

// return A-B, where A is *this
// this operator invokes elementwise_addition()
Matrix Matrix::operator-(const Matrix &B) const {
    // exception handling
    if (nrows!=B.Nrows() || ncols!=B.Ncols()) {
        *Basics::err << "Matrix::operator-(const Matrix&) const: matrix dimensions must agree!" << endl;
        Basics::quit(1);
    }

    // pointer for output
    Real *sc = NULL;

    // compute C = A - B, where A is *this and C is stored in sc[]
    elementwise_addition(nrows*ncols, Store(), -1.0, B.Store(), sc);

    // return the result
    return Matrix(nrows, ncols, sc);
}

// return a const reference to A after computing A = A+B, where A is *this
// this operator invokes elementwise_addition()
const Matrix& Matrix::operator+=(const Matrix &B) {
    // exception handling
    if (nrows!=B.Nrows() || ncols!=B.Ncols()) {
        *Basics::err << "Matrix::operator+=(const Matrix&): matrix dimensions must agree!" << endl;
        Basics::quit(1);
    }

    // compute A += B, where A is *this
    elementwise_addition(nrows*ncols, Store(), 1.0, B.Store(), store);

    // return the result
    return *this;
}

// return a const reference to A after computing A = A-B, where A is *this
// this operator invokes elementwise_addition()
const Matrix& Matrix::operator-=(const Matrix &B) {
    // exception handling
    if (nrows!=B.Nrows() || ncols!=B.Ncols()) {
        *Basics::err << "Matrix::operator-=(const Matrix&): matrix dimensions must agree!" << endl;
        Basics::quit(1);
    }

    // compute A -= B, where A is *this
    elementwise_addition(nrows*ncols, Store(), -1.0, B.Store(), store);

    // return the result
    return *this;
}


// matrix - symmetric matrix addition/subtraction
// return A+B, where A is *this
// this operator invokes matrix_symmatrix_addition()
Matrix Matrix::operator+(const SymmetricMatrix &B) const {
    // exception handling
    if (nrows!=ncols || nrows!=B.Nrows()) {
        // note that B.Nrows()==B.Ncols() since B is symmetric
        *Basics::err << "Matrix::operator+(const SymmetricMatrix&) const: matrix dimensions must agree!" << endl;
        Basics::quit(1);
    }

    // pointer for output
    Real *sc = NULL;

    // compute (*this)+B
    matrix_symmatrix_addition(nrows, Store(), 1.0, B.Store(), sc);

    // return the result
    return Matrix(nrows, nrows, sc);
}

// return A-B, where A is *this
// this operator invokes matrix_symmatrix_addition()
Matrix Matrix::operator-(const SymmetricMatrix &B) const {
    // exception handling
    if (nrows!=ncols || nrows!=B.Ncols()) {
        // note that B.Nrows()==B.Ncols() since B is symmetric
        *Basics::err << "Matrix::operator-(const SymmetricMatrix&) const: matrix dimensions must agree!" << endl;
        Basics::quit(1);
    }

    // pointer for output
    Real *sc = NULL;

    // compute (*this)-B
    matrix_symmatrix_addition(nrows, Store(), -1.0, B.Store(), sc);

    // return the result
    return Matrix(nrows, ncols, sc);
}

// return a const reference to A after computing A=A+B, where A is *this
// this operator invokes matrix_symmatrix_addition()
const Matrix& Matrix::operator+=(const SymmetricMatrix &B) {
    // exception handling
    if (nrows!=B.Nrows() || ncols!=B.Ncols()) {
        *Basics::err << "Matrix::operator+=(const SymmetricMatrix&) const: matrix dimensions must agree!" << endl;
        Basics::quit(1);
    }

    // compute *this = *this + B
    matrix_symmatrix_addition(nrows, Store(), 1.0, B.Store(), store);

    // return the result
    return *this;
}

// return a const reference to A after computing A=A-B, where A is *this
// this operator invokes matrix_symmatrix_addition()
const Matrix& Matrix::operator-=(const SymmetricMatrix &B) {
    // exception handling
    if (nrows!=B.Nrows() || ncols!=B.Ncols()) {
        *Basics::err << "Matrix::operator-=(const SymmetricMatrix&) const: matrix dimensions must agree!" << endl;
        Basics::quit(1);
    }

    // compute *this = *this - B
    matrix_symmatrix_addition(nrows, Store(), -1.0, B.Store(), store);

    // return the result
    return *this;
}


// matrix - sparse matrix addition/subtraction
// return A+B, where A is *this
// this operator invokes matrix_spmatrix_addition()
Matrix Matrix::operator+(const SparseMatrix &B) const {
    // exception handling
    if (nrows != B.Nrows() || ncols != B.Ncols()) {
        *Basics::err << "Matrix::operator+(const SparseMatrix&) const: matrix dimensions must agree!" << endl;
        Basics::quit(1);
    }
    // pointer for output
    Real *sc = NULL;

    // compute C = A+B, where A is *this and C is stored in sc[]
    matrix_spmatrix_addition(nrows, ncols, Store(), 1.0, B.Store(), B.RowIndex(), B.ColPointer(), sc);

    // return the result
    return Matrix(nrows, ncols, sc);
}

// return A-B, where A is *this
// this operator invokes matrix_spmatrix_addition()
Matrix Matrix::operator-(const SparseMatrix &B) const {
    // exception handling
    if (nrows != B.Nrows() || ncols != B.Ncols()) {
        *Basics::err << "Matrix::operator-(const SparseMatrix&) const: matrix dimensions must agree!" << endl;
        Basics::quit(1);
    }
    // pointer for output
    Real *sc = NULL;

    // compute C = A-B, where A is *this and C is stored in sc[]
    matrix_spmatrix_addition(nrows, ncols, Store(), -1.0, B.Store(), B.RowIndex(), B.ColPointer(), sc);

    // return the result
    return Matrix(nrows, ncols, sc);
}

// return a const reference to A after computing A=A+B, where A is *this
// this operator invokes matrix_spmatrix_addition()
const Matrix& Matrix::operator+=(const SparseMatrix &B) {
    // exception handling
    if (nrows != B.Nrows() || ncols != B.Ncols()) {
        *Basics::err << "Matrix::operator+=(const SparseMatrix&): matrix dimensions must agree!" << endl;
        Basics::quit(1);
    }
    // compute A += B, where A is *this
    matrix_spmatrix_addition(nrows, ncols, Store(), 1.0, B.Store(), B.RowIndex(), B.ColPointer(), store);

    // return the result
    return *this;
}

// return a const reference to A after computing A=A-B, where A is *this
// this operator invokes matrix_spmatrix_addition()
const Matrix& Matrix::operator-=(const SparseMatrix &B) {
    if (nrows != B.Nrows() || ncols != B.Ncols()) {
        *Basics::err << "Matrix::operator-=(const SparseMatrix&): matrix dimensions must agree!" << endl;
        Basics::quit(1);
    }
    // compute A -= B, where A is *this
    matrix_spmatrix_addition(nrows, ncols, Store(), -1.0, B.Store(), B.RowIndex(), B.ColPointer(), store);

    // return a const reference to *this
    return *this;
}


// matrix-vector multiplication
// return A*v, where A is *this
// this operator invokes matrix_vector_multiplication()
Vector Matrix::operator*(const Vector &v) const {
    if (ncols != v.Length()) {
        *Basics::err << "Matrix::operator*(const Vector&) const: inner matrix-vector dimensions must agree!" << endl;
        Basics::quit(1);
    }

#ifdef USE_MWBLAS
    // MWBLAS routines sgemv and dgemv require non-const storage for the input matrix and vector,
    // so it cannot be incorporated into matrix_vector_multiplication and therefore is put here
    Real *sw = new Real[nrows];
    char trans = 'N';
    Real alp = 1.0, bet = 0.0;
    mwSignedIndex inc = 1;
    mwSignedIndex nrows2 = (mwSignedIndex)nrows;
    mwSignedIndex ncols2 = (mwSignedIndex)ncols;
    #ifdef USE_SINGLE
        sgemv(&trans, &nrows2, &ncols2, &alp, store, &nrows2, v.store, &inc, &bet, sw, &inc);
    #else
        dgemv(&trans, &nrows2, &ncols2, &alp, store, &nrows2, v.store, &inc, &bet, sw, &inc);
    #endif
#else
    // pointer for output
    Real *sw = NULL;

    // compute (*this)*v
    matrix_vector_multiplication(nrows, ncols, Store(), v.Store(), sw);
#endif

    // return the result
    return Vector(nrows, sw);
}

// matrix - matrix multiplication
// return A*B, where A is *this
// this operator invokes matrix_matrix_multiplication()
Matrix Matrix::operator*(const Matrix &B) const {
    if (ncols != B.Nrows()) {
        *Basics::err << "Matrix::operator*(const Matrix&) const: inner matrix dimensions must agree!" << endl;
        Basics::quit(1);
    }
#ifdef USE_MWBLAS
    // MWBLAS routines sgemm and dgemm require non-const storage for reading input matrices,
    // so it cannot be incorporated into matrix_matrix_multiplication and therefore is put here
    Real *sc = new Real[nrows*B.Ncols()];
    char trans = 'N';
    Real alp = 1.0, bet = 0.0;
    mwSignedIndex m = (mwSignedIndex)nrows;
    mwSignedIndex k = (mwSignedIndex)ncols;
    mwSignedIndex n = (mwSignedIndex)B.Ncols();
    #ifdef USE_SINGLE
        sgemm(&trans, &trans, &m, &n, &k, &alp, store, &m, B.store, &k, &bet, sc, &m);
    #else
        dgemm(&trans, &trans, &m, &n, &k, &alp, store, &m, B.store, &k, &bet, sc, &m);
    #endif
#else
    // pointer for output
    Real *sc = NULL;

    // compute C = A*B, where A is *this, and C is stored in sc[]
    matrix_matrix_multiplication(nrows, ncols, B.Ncols(), Store(), B.Store(), sc);
#endif

    // return the result
    return Matrix(nrows, B.Ncols(), sc);
}

// return a const reference to A after computing A = A*B, where A is *this
// this operator (eventually) invokes matrix_matrix_multiplication()
const Matrix& Matrix::operator*=(const Matrix &B) {
    // exception handling
    if (ncols != B.Nrows()) {
        *Basics::err << "Matrix::operator*=(const Matrix&): inner matrix dimensions must agree!" << endl;
        Basics::quit(1);
    }

    // compute A = A*B and return the result, where A is this
    return *this = (*this)*B;
}


// matrix - symmetric matrix multiplication
// return A*B, where A is *this
// this operator invokes matrix_symmatrix_multiplication()
Matrix Matrix::operator*(const SymmetricMatrix &B) const {
    // exception handling
    if (ncols != B.Nrows()) {
        *Basics::err << "Matrix::operator*(const SymmetricMatrix&) const: inner matrix dimensions must agree!" << endl;
        Basics::quit(1);
    }

    // pointer for output
    Real *sc = NULL;

    // compute C = A*B, where A is *this, and C is stored in sc[]
    matrix_symmatrix_multiplication(nrows, ncols, Store(), B.Store(), sc);

    // return the result
    return Matrix(nrows, ncols, sc);
}

// return a const reference to A after computing A = A*B, where A is *this
// this operator (eventually) invokes matrix_symmatrix_multiplication()
const Matrix& Matrix::operator*=(const SymmetricMatrix &B) {
    // exception handling
    if (ncols != B.Nrows()) {
        *Basics::err << "Matrix::operator*=(const SymmetricMatrix&): inner matrix dimensions must agree!" << endl;
        Basics::quit(1);
    }

    // return a const reference to A after computing A = A*B, where A is *this
    return *this = (*this)*B;
}


// matrix - sparse matrix multiplication
// return A*B, where A is *this
// this operator invokes matrix_symmatrix_multiplication()
Matrix Matrix::operator*(const SparseMatrix &B) const {
    if (ncols != B.Nrows()) {
        *Basics::err << "Matrix::operator*(const SparseMatrix&) const: inner matrix dimensions must agree!" << endl;
        Basics::quit(1);
    }
    // pointer for output
    Real *sc = NULL;

    // compute C = A*B, where A is *this and C is stored in sc[]
    matrix_spmatrix_multiplication(nrows, ncols, B.Ncols(), Store(), B.Store(), B.RowIndex(), B.ColPointer(), sc);

    // return the result
    return Matrix(nrows, B.Ncols(), sc);
}

// return a const reference to A after computing A = A*B, where A is *this
// this operator (eventually) invokes matrix_symmatrix_multiplication()
const Matrix& Matrix::operator*=(const SparseMatrix &B) {
    // exception handling
    if (ncols != B.Nrows()) {
        *Basics::err << "Matrix::operator*=(const SparseMatrix&): inner matrix dimensions must agree!" << endl;
        Basics::quit(1);
    }

    // return a reference to A after computing A = A*B, where A is *this
    return *this = (*this)*B;
}


// permute columns of A, such that
//    B(i,perm[j-1])   = A(i,j) for i=1,...,nrows and j=1,...,ncols (if indexFrom1==true), or
//    B(i,perm[j-1]+1) = A(i,j) for i=1,...,nrows and j=1,...,ncols (if indexFrom1==false),
// where A is *this and B is the return matrix
// the indices perm[0,...,ncols-1] range over 0,1,...,ncols-1 (if indexFrom1==false) or over 1,2,...,ncols (if indexFrom1==true)
// this member function invokes permute_matrix_columns()
// this member function resembles B(:,perm) = A in OCTAVE/MATLAB
Matrix Matrix::permuteColumns(const mkIndex *perm, bool indexFrom1) const {
    // pointer for output
    Real *sb = NULL;

    // compute the permutation
    permute_matrix_columns(nrows, ncols, perm, Store(), sb, indexFrom1);

    // return the result
    return Matrix(nrows, ncols, sb);
}

// permute rows of A, such that
//    B(perm[i-1],j)   = A(i,j) for i=1,...,nrows and j=1,...,ncols (if indexFrom1==true), or
//    B(perm[i-1]+1,j) = A(i,j) for i=1,...,nrows and j=1,...,ncols (if indexFrom1==false),
// where A is *this and B is the return matrix
// the indices perm[0,...,nrows-1] range over 0,1,...,nrows-1 (if indexFrom1==false) or over 1,2,...,nrows (if indexFrom1==true)
// this member function invokes permute_matrix_rows()
// this member function resembles B(perm,:) = A in OCTAVE/MATLAB
Matrix Matrix::permuteRows(const mkIndex *perm, bool indexFrom1) const {
    // pointer for output
    Real *sb = NULL;

    // compute the permutation
    permute_matrix_rows(nrows, ncols, perm, Store(), sb, indexFrom1);

    // return the result
    return Matrix(nrows, ncols, sb);
}

// permute rows and columns of A, such that
//    B(rperm[i-1],cperm[j-1])     = A(i,j) for i=1,...,nrows and j=1,...,ncols (if indexFrom1==true), or
//    B(rperm[i-1]+1,cperm[j-1]+1) = A(i,j) for i=1,...,nrows and j=1,...,ncols (if indexFrom1==false),
// where A is *this and B is the return matrix
// if cperm==NULL, then symmetric permutation will be applied (i.e., rperm[] will also be used for cperm[])
// the indices rperm[0,...,nrows-1] range over 0,1,...,nrows-1 (if indexFrom1==false) or over 1,2,...,nrows (if indexFrom1==true)
// the indices cperm[0,...,ncols-1] range over 0,1,...,ncols-1 (if indexFrom1==false) or over 1,2,...,ncols (if indexFrom1==true)
// this member function (eventually) invokes permute_matrix_rows() and permute_matrix_columns()
// this member function resembles B(rperm,cperm) = A in OCTAVE/MATLAB
Matrix Matrix::permuteRowsAndColumns(const mkIndex *rperm, const mkIndex *cperm, bool indexFrom1) const {
    if (cperm == NULL) {
        // exception handling
        if (nrows != ncols) {
            *Basics::err << "Matrix::permuteRowsAndColumns(const mkIndex*, const mkIndex*, bool) const: symmetric permutation requires the matrix being square!" << endl;
            Basics::quit(1);
        }
        // symmetrically permute rows and columns
        return permuteRows(rperm, indexFrom1).permuteColumns(rperm, indexFrom1);
    }
    // permute rows and columns
    return permuteRows(rperm, indexFrom1).permuteColumns(cperm, indexFrom1);
}


#ifdef USE_MEX
// convert *mxArray to Matrix
Matrix Matrix::mxArray2Matrix(const mxArray *mA) {
    // number of rows and number of columns
    mkIndex nr, nc;
    // pointer for output
    Real *sa = NULL;

    // convert
    mxArray_to_matrix(mA, nr, nc, sa);

    // return the result
    return Matrix(nr, nc, sa);
}

// convert Matrix to *mxArray
mxArray *Matrix::Matrix2mxArray(const Matrix &A) {
    return matrix_to_mxArray(A.Nrows(), A.Ncols(), A.Store());
}
#endif


#ifdef USE_NAMESPACE
}  // end of namespace NEWMAT
#endif
