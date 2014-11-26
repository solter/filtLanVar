//============================================================================
// the MATKIT class SparseMatrix
// coded by H.-r. Fang
// last update November, 2011
//============================================================================

#include <stdio.h>   // for fopen, fscanf, fclose
#include <stdlib.h>  // for exit, atof, atoi, atol, srand, rand, etc.
#include <string.h>  // for memcpy, strlen, strcmp, strncmp, etc.
#include <stddef.h>  // for NULL pointer, pointer subtraction, etc.
#include <ctype.h>   // for tolower
#include <math.h>    // for sqrt, pow, etc.
#include <iostream>  // for cout, cerr, endl, etc. (under namespace std)

#include "matkit.h"

using std::endl;

#ifdef USE_NAMESPACE
namespace MATKIT {
#endif


Real SparseMatrix::memoryExpansionFactor = 1.2;


// display sparse matrix A with operator <<
// in other words, print out sparse matrix A to ostream s
// this routine invokes print_spmatrix()
std::ostream& operator<<(std::ostream &s, const SparseMatrix &A) {
    const mkIndex *colPtr = A.ColPointer();
    const mkIndex *rowIdx = A.RowIndex();
    const Real *str = A.Store();
    mkIndex shft = colPtr[0];
    mkIndex j = 0;
    std::ostream *sptr = &s;

    #ifdef USE_MEX
    if (&s == &std::cout || &s == &std::cerr) {
        // stdout and stderr are not supported by the mex interface
        // use mexostream (which inherits ostream) "mexcout" for printing messages to OCTAVE/MATLAB console
        // see matkitdef.h for more information
        sptr = &mexcout;
    }
    #endif
    *sptr << "compressed sparse columns (rows = " << A.Nrows() << ", columns = " << A.Ncols() << ", nnz = " << colPtr[A.Ncols()]-colPtr[0] << ")" << endl;
    for (mkIndex i=shft; i<colPtr[A.Ncols()]; i++) {
        while (colPtr[j] <= i)
            j ++;
        *sptr << "    (" << rowIdx[i-shft]+1u << "," << j << ")";
        *sptr << "     " << str[i-shft] << endl;
    }

    return s;
}

// read a sparse matrix from a matrix-market file and return the result
// fileName[] is the file name of the matrix-market file
// sort       tells whether to sort elements in each column with respect to the row indices
// usually a matrix-market file has elements sorted, in which case sort==true is unnecessary
// this member function invokes matrix_market_spmatrix_read()
SparseMatrix SparseMatrix::mmread(const char fileName[], bool sort) {
    // number of rows and number of columns
    mkIndex nrows, ncols;

    // variables for coordinate format
    mkIndex nnz;
    mkIndex *ridx = NULL, *cidx = NULL;
    Real *sa = NULL;
    char symm;

    // read the matrix in coordinate format
    int info = matrix_market_spmatrix_read(fileName, nrows, ncols, nnz, sa, ridx, cidx, symm);

    // in case of any exception, report and quit
    switch (info) {
        case -1:  // file open error
            *Basics::err << "SparseMatrix::mmread(const char*, bool): file open error!" << endl;
            Basics::quit(1);
        case -2:  // file is empty
            *Basics::err << "SparseMatrix::mmread(const char*, bool): file is empty!" << endl;
            Basics::quit(1);
        case -3:  // the first line in the input file is empty
            *Basics::err << "SparseMatrix::mmread(const char*, bool): the first line of the input file is empty!" << endl;
            Basics::quit(1);
        case -4:  // invalid header
            *Basics::err << "SparseMatrix::mmread(const char*, bool): invalid header!" << endl;
            Basics::quit(1);
        case -5:  // invalid number of rows, number of columns, or number of nonzero elements
            *Basics::err << "SparseMatrix::mmread(const char*, bool): invalid number of rows, number of columns, or number of nonzero elements!" << endl;
            Basics::quit(1);
        case -6:  // matrix data is incomplete or invalid
            *Basics::err << "SparseMatrix::mmread(const char*, bool): the matrix data is incomplete or invalid!" << endl;
            Basics::quit(1);
        case  1:  // the input file is a dense matrix
            *Basics::err << "SparseMatrix::mmread(const char*, bool): this routine does not read a dense matrix!" << endl;
            Basics::quit(1);
        case  2:  // the input file is a complex (sparse) matrix
            *Basics::err << "SparseMatrix::mmread(const char*, bool): this routine does not read a complex (sparse) matrix!" << endl;
            Basics::quit(1);
        case  3:  // the intput file contains only sparsity information
            *Basics::err << "SparseMatrix::mmread(const char*, bool): the input file contains only sparsity information!" << endl;
            Basics::quit(1);
        default:  // normal case
            if (info != 0) {
                *Basics::err << "SparseMatrix::mmread(const char*, bool): invalid return flag of matrix_market_spmatrix_read()!" << endl;
                Basics::quit(1);
            }
    }

    // pointers for output in CSC format
    Real *store = NULL;
    mkIndex *rowIdx = NULL;
    mkIndex *colPtr = NULL;

    // convert coordinates to CSC format
    if (coordinates_to_csc(ncols, nnz, sa, ridx, cidx, store, rowIdx, colPtr, symm) != 0) {
        // shouldn't happen; exceptions, if any, should have been handled with the return flag of matrix_market_spmatrix_read()
        *Basics::err << "SparseMatrix::mmread(const char*, bool): conversion from coordinate format to CSC format by coordinates_to_csc() failed!" << endl;
        Basics::quit(1);
    }

    // free memory for the matrix in coordinate format
    delete [] sa;
    delete [] ridx;
    delete [] cidx;

    // sort the elements in each column with respect to the row indices
    if (sort)
        sort_csc_elements_in_each_column(nrows, ncols, store, rowIdx, colPtr);

    // return the result
    return SparseMatrix(nrows, ncols, store, rowIdx, colPtr);
}


// return an n-by-n (symmetric) matrix with ones on the main diagonal and zeros elsewhere
// this static function resembles speye in OCTAVE/MATLAB
SparseMatrix SparseMatrix::eye(mkIndex nr, mkIndex nc) {
    // allocate memory
    mkIndex mrc = (nr<=nc ? nr : nc);
    Real *store = new Real[mrc];
    mkIndex *rowIdx = new mkIndex[mrc];
    mkIndex *colPtr = new mkIndex[nc+1u];

    // set the sparse eye
    for (mkIndex i=0; i<mrc; i++) {
        store[i] = 1.0;
        rowIdx[i] = i;
        colPtr[i] = i;
    }
    for (mkIndex i=mrc; i<=nc; i++)
        colPtr[i] = mrc;

    // return the result
    return SparseMatrix(nr, nc, store, rowIdx, colPtr);
}

// return a tridiagonal matrix A as a sparse matrix with three diagonals ldiag[], diag[], and udiag[]
// the lower subdiagonal, the main diagonal, and the upper subdiagonal of A are formed by
// ldiag[0,...,nrc-2], diag[0,...,nrc-1], and udiag[0,...,nrc-2], respectively
// to be specific, A(i+1,i) = ldiag[i-1] for i=1,...,nrc-1,
//                 A(i,i  ) =  diag[i-1] for i=1,...,nrc,
//             and A(i,i+1) = udiag[i]   for i=1,...,nrc-1
// this static function invokes tridiagonal_to_spmatrix()
SparseMatrix SparseMatrix::tridiagonalMatrix(mkIndex nrc, const Real *ldiag, const Real *mdiag, const Real *udiag) {
    // special case, an empty matrix
    if (nrc == 0)
        return SparseMatrix();

    // pointers for output
    Real *sa = NULL;
    mkIndex *rowIdx = NULL;
    mkIndex *colPtr = NULL;

    // form the sparse matrix
    tridiagonal_to_spmatrix(nrc, ldiag, mdiag, udiag, sa, rowIdx, colPtr);

    // return the result
    return SparseMatrix(nrc, nrc, sa, rowIdx, colPtr);
}


// constructors
// a constructor for an empty matrix
SparseMatrix::SparseMatrix() {
    nrows = 0;
    ncols = 0;
    store = NULL;
    rowIndex = NULL;
    colPointer = NULL;
    maxnnz = 0;
}

// a constructor for an nr-by-nc (sparse) matrix of zeros
SparseMatrix::SparseMatrix(mkIndex nr, mkIndex nc) {
    nrows = nr;
    ncols = nc;
    store = NULL;
    rowIndex = NULL;
    colPointer = new mkIndex[nc+1u];
    for (mkIndex j=0; j<=nc; j++)
        colPointer[j] = 0;
    maxnnz = 0;
}

// a constructor for an nr-by-nc (sparse) matrix stored in CSC format in store0[], rowIndex0[], colPointer0[]
// note that it performs a shallow copy of store0, rowIndex0, and colPointer0, and the memory will be freed by the destructor
// if maxnnz0 is greater than the number of nonzero elements, then it is assumed that
// memory of size maxnnz0 is allocated, with the address pointed to by store0 and so is it for rowIndex0
SparseMatrix::SparseMatrix(mkIndex nr, mkIndex nc, Real *store0, mkIndex *rowIndex0, mkIndex *colPointer0, mkIndex maxnnz0) {
    // set dimensions
    nrows = nr;
    ncols = nc;

    store = store0;
    rowIndex = rowIndex0;
    colPointer = colPointer0;

    maxnnz = colPointer[ncols]-colPointer[0];
    if (maxnnz < maxnnz0)
        maxnnz = maxnnz0;
    // the default maxnnz0 is 0, in which case we have maxnnz = colPointer[ncols]-colPointer[0]
}

// a copy constructor
// this constructor invokes memory_xcopy() for copying nonzero elements
SparseMatrix::SparseMatrix(const SparseMatrix &B) {
    // set dimensions
    nrows = B.Nrows();
    ncols = B.Ncols();
    maxnnz = B.colPointer[ncols]-colPointer[0];

    // allocate memory and copy the elements
    store = new Real[maxnnz];
    memory_xcopy(maxnnz, B.Store(), store);

    // allocate memory and copy the row indices
    rowIndex = new mkIndex[maxnnz];
    memcpy(rowIndex, B.RowIndex(), maxnnz*sizeof(mkIndex));

    // allocate memory and copy the column pointers
    colPointer = new mkIndex[ncols+1u];
    memcpy(colPointer, B.ColPointer(), (ncols+1u)*sizeof(mkIndex));
}

// a destructor
SparseMatrix::~SparseMatrix() {
    delete [] store;
    delete [] rowIndex;
    delete [] colPointer;
}


// convert *this to a general matrix
// this static function invokes sparse_to_full()
// this static function resembles full(A) in OCTAVE/MATLAB, where A is *this
Matrix SparseMatrix::toFull() const {
    // pointer for output
    Real *sb = NULL;

    // convert
    sparse_to_full(nrows, ncols, Store(), RowIndex(), ColPointer(), sb);

    // return the result
    return Matrix(nrows, ncols, sb);
}

// set *this as a sparse matrix of size nr-by-nc stored in CSC format in store0[], rowIndex0[], colPointer0[]
// note that it performs a shallow copy of store0, rowIndex0, and colPointer0, and the memory will be freed by the destructor
// if maxnnz0 is greater than the number of nonzero elements, then it is assumed that
// memory of size maxnnz0 is allocated, with the address pointed to by store0 and so is it for rowIndex0
void SparseMatrix::set(mkIndex nr, mkIndex nc, Real *store0, mkIndex *rowIndex0, mkIndex *colPointer0, mkIndex maxnnz0) {
    // set dimensions
    nrows = nr;
    ncols = nc;
    if (store != store0) {
        // free memory and do a shallow copy
        delete [] store;
        store = store0;
    }
    if (rowIndex != rowIndex0) {
        // free memory and do a shallow copy
        delete [] rowIndex;
        rowIndex = rowIndex0;
    }
    if (colPointer != colPointer0) {
        // free memory and do a shallow copy
        delete [] colPointer;
        colPointer = colPointer0;
    }

    maxnnz = colPointer[ncols]-colPointer[0];
    if (maxnnz < maxnnz0)
        maxnnz = maxnnz0;
}


// a submatrix formed by columns j1,...,j2 of *this
// this member function invokes submatrix_of_spmatrix()
// this member function resembles A(:,j1:j2) in OCTAVE/MATLAB, where A is *this
SparseMatrix SparseMatrix::columns(mkIndex j1, mkIndex j2) const {
    // pointers for the output submatrix
    Real *store2 = NULL;
    mkIndex *rowIdx2 = NULL;
    mkIndex *colPtr2 = NULL;

    // compute the submatrix
    submatrix_of_spmatrix(ncols, store, rowIndex, colPointer, 1, nrows, j1, j2, store2, rowIdx2, colPtr2);

    // number of columns of the submatrix
    mkIndex nc2 = (j1<=j2) ? (j2-j1+1u) : (j1-j2+1u);

    // return the result
    return SparseMatrix(nrows, nc2, store2, rowIdx2, colPtr2);
}

// a submatrix formed by rows i1,...,i2 of *this
// this member function invokes submatrix_of_spmatrix()
// this member function resembles A(i1:i2,:) in OCTAVE/MATLAB, where A is *this
SparseMatrix SparseMatrix::rows(mkIndex i1, mkIndex i2) const {
    // pointers for the output submatrix
    Real *store2 = NULL;
    mkIndex *rowIdx2 = NULL;
    mkIndex *colPtr2 = NULL;

    // compute the submatrix
    submatrix_of_spmatrix(ncols, store, rowIndex, colPointer, i1, i2, 1, ncols, store2, rowIdx2, colPtr2);

    // number of rows of the submatrix
    mkIndex nr2 = (i1<=i2) ? (i2-i1+1u) : (i1-i2+1u);

    // return the result
    return SparseMatrix(nr2, ncols, store2, rowIdx2, colPtr2);
}

// a submatrix formed by elements in rows i1,...,i2 and columns j1,...,j2 of *this
// this member function invokes submatrix_of_matrix()
// this member function resembles A(i1:i2,j1:j2) in OCTAVE/MATLAB, where A is *this
SparseMatrix SparseMatrix::subMatrix(mkIndex i1, mkIndex i2, mkIndex j1, mkIndex j2) const {
    // pointers for the output submatrix
    Real *store2 = NULL;
    mkIndex *rowIdx2 = NULL;
    mkIndex *colPtr2 = NULL;

    // compute the submatrix
    submatrix_of_spmatrix(ncols, store, rowIndex, colPointer, i1, i2, j1, j2, store2, rowIdx2, colPtr2);

    // numbers of rows and columns of the submatrix
    mkIndex nr2 = (i1<=i2) ? (i2-i1+1u) : (i1-i2+1u);
    mkIndex nc2 = (j1<=j2) ? (j2-j1+1u) : (j1-j2+1u);

    // return the result
    return SparseMatrix(nr2, nc2, store2, rowIdx2, colPtr2);
}


// a vector formed by the elements in the k-th diagonal of *this in order
// this member function invokes spmatrix_diag()
// this member function resembles diag(A,k) in OCTAVE/MATLAB
Vector SparseMatrix::diag(mkSignedIndex k) const {
    // pointer for output
    Real *sv = NULL;

    // extract the k-th diagonal
    mkIndex len = spmatrix_diag(nrows, ncols, store, rowIndex, colPointer, sv, k);

    // return the result
    return Vector(len, sv);
}

// sparse matrix transpose
// this member function invokes spmatrix_transpose()
// this member function resembles A' or equivalently transpose(A) in OCTAVE/MATLAB, where A is *this
SparseMatrix SparseMatrix::transpose() const {
    Real *store2 = NULL;
    mkIndex *rowIdx2 = NULL;
    mkIndex *colPtr2 = NULL;
    spmatrix_transpose(nrows, ncols, store, rowIndex, colPointer, store2, rowIdx2, colPtr2);

    // return the transpose of *this
    return SparseMatrix(ncols, nrows, store2, rowIdx2, colPtr2);
}


// Frobenius norm of *this
Real SparseMatrix::normFrobenius() const {
#if defined(USE_MWBLAS)
    mwSignedIndex nnz = (mwSignedIndex)(colPointer[ncols]-colPointer[0]);
    mwSignedIndex inc = 1;
    #ifdef USE_SINGLE
        return snrm2(&nnz, store, &inc);
    #else
        return dnrm2(&nnz, store, &inc);
    #endif
#elif defined(USE_BLAS)
    int nnz = (int)(colPointer[ncols]-colPointer[0]);
    int inc = 1;
    #ifdef USE_SINGLE
        return snrm2_(&nnz, store, &inc);
    #else
        return dnrm2_(&nnz, store, &inc);
    #endif
#else
    return sqrt(sumSquare());
#endif
}

// sum of squared elements of *this
Real SparseMatrix::sumSquare() const {
    Real sum = 0.0;
    Real *s = store;
    mkIndex nnz = colPointer[ncols]-colPointer[0];
    while (nnz--)
        sum += Basics::square(*s++);
    return sum;
}

// 1-norm of *this
Real SparseMatrix::norm1() const {
    Real nrm = 0.0;
    Real val = 0.0;
    Real *sa = store;
    mkIndex j = 0;
    mkIndex z = colPointer[ncols];
    for (mkIndex i=colPointer[0]; i<z; i++) {
        while (colPointer[j] <= i) {
            // at the end of the current column, check whether the current column 1-norm is greater than previous ones
            // if it is, update nrm with the current column 1-norm
            nrm = (val>nrm ? val : nrm);
            val = 0.0;
            j ++;
        }
        val += Basics::abs(*sa++);
    }
    return nrm;
}

// infinity-norm of *this
Real SparseMatrix::normInfinity() const {
    // allocate work space of size nrows and initialize it as zero
    Real *work = new Real[nrows];
    Real *w = work;
    mkIndex ii = nrows;
    while (ii--)
        *w++ = 0.0;

    // compute the infinity-norm of each row
    Real *sa = store;
    mkIndex *ri = rowIndex;
    mkIndex jj = colPointer[ncols]-colPointer[0];
    while (jj--)
        work[*ri++] += Basics::abs(*sa++);

    // the matrix infinity-norm is the maximum row infinity-norm
    Real nrm = 0.0;
    for (mkIndex i=0; i<nrows; i++)
        nrm = (nrm>=work[i] ? nrm : work[i]);

    // free the work space
    delete [] work;

    // return the result
    return nrm;
}


// static function setMemoryExpansionFactor
// the default memoryExpansionFactor is 1.2, shown in spmatrix.cpp
// the memoryExpansionFactor is used by SparseMatrix::operator(mkIndex,mkIndex)
void SparseMatrix::setMemoryExpansionFactor(Real newExpansionFactor) {
    if (newExpansionFactor <= 1.0 || newExpansionFactor > 2.0) {
        *Basics::err << "SparseMatrix::setMemoryExpansionFactor(Real): expansion factor should be >1.0 and <= 2.0!" << endl;
        Basics::quit(1);
    }

    // set the memory expansion factor
    memoryExpansionFactor = newExpansionFactor;
}

// access A(i,j)
// return a reference to A(i,j), i.e. can be used to write A(i,j), where A is *this
// if A(i,j)==0, then a space will be created for it
Real& SparseMatrix::operator()(mkIndex i, mkIndex j) {
    // exception handling
    if (i==0 || j==0) {
        *Basics::err << "SparseMatrix::operator()(mkIndex, mkIndex): subscript indices must be positive!" << endl;
        Basics::quit(1);
    }
    if (i>nrows || j>ncols) {
        *Basics::err << "SparseMatrix::operator()(mkIndex, mkIndex): index exceeds matrix dimensions!" << endl;
        Basics::quit(1);
    }

    mkIndex i1 = i-1;
    mkIndex k;
    for (k=colPointer[j-1u]-colPointer[0]; k<colPointer[j]-colPointer[0]; k++) {
        if (rowIndex[k] == i1) {
            // A(i,j) is nonzero (or occupies a space at least); return the reference to it
            return store[k];
        }
        if (rowIndex[k] > i1) {
            // A(i,j) is zero; have to create a space for it
            break;
        }
    }
    mkIndex nnz = colPointer[ncols]-colPointer[0];
    if (maxnnz == nnz) {
        // no more space; expand
        mkIndex old_maxnnz = maxnnz;
        maxnnz *= memoryExpansionFactor;
        if (maxnnz <= old_maxnnz)
            maxnnz = maxnnz + 1u;
        Real *store2 = new Real[maxnnz];
        mkIndex *rowIndex2 = new mkIndex[maxnnz];
        // copy elements in 2 installments
        if (k > 0) {
            memory_xcopy(k, Store(), store2);
            memcpy(rowIndex2, RowIndex(), k*sizeof(mkIndex));
        }
        if (k < nnz) {
            memory_xcopy(nnz-k, Store()+k, store2+k+1);
            memcpy(rowIndex2+k+1, RowIndex()+k, (nnz-k)*sizeof(mkIndex));
        }
        // free memory and reset the pointers
        delete [] store;
        delete [] rowIndex;
        store = store2;
        rowIndex = rowIndex2;
    }
    else {
        // shift the elements / row indices
        mkIndex jj = nnz-k;
        Real *s2 = store+nnz;
        Real *s1 = s2-1;
        mkIndex *r2 = rowIndex-k;
        mkIndex *r1 = r2-1;
        while (jj--) {
            *s2-- = *s1--;
            *r2-- = *r1--;
        }
        // after adding the new entry, nnz is increased by 1
    }
    // increase the value of column pointer j0 by 1 for j0=j,...,ncols
    for (mkIndex j0=j; j0<=ncols; j0++)
        colPointer[j0] ++;

    // the entries of A(i,j)
    store[k] = 0.0;
    rowIndex[k] = i1;
    // return a reference to A(i,j)
    return store[k];
}

// read A(i,j) via a const reference or pointer to A, where A is *this
Real SparseMatrix::operator()(mkIndex i, mkIndex j) const {
    // exception handling
    if (i==0 || j==0) {
        *Basics::err << "SparseMatrix::operator()(mkIndex, mkIndex) const: subscript indices must be positive!" << endl;
        Basics::quit(1);
    }
    if (i>nrows || j>ncols) {
        *Basics::err << "SparseMatrix::operator()(mkIndex, mkIndex) const: index exceeds matrix dimensions!" << endl;
        Basics::quit(1);
    }

    // return the value of A(i,j)
    mkIndex i1 = i-1;
    for (mkIndex k=colPointer[j-1u]-colPointer[0]; k<colPointer[j]-colPointer[0]; k++) {
        if (rowIndex[k] == i1)
            return store[k];
        if (rowIndex[k] > i1)
            break;
    }
    return 0.0;
}


// operator =
// return a const reference to A after setting A=B, where A is *this
// this operator invokes memory_xcopy() for copying nonzero elements
const SparseMatrix& SparseMatrix::operator=(const SparseMatrix &B) {
    // number of nonzero elements
    mkIndex nnz = B.colPointer[B.ncols]-B.colPointer[0];

    if (maxnnz < nnz) {
        // memory is inefficient, so re-allocate
        maxnnz = nnz;
        delete [] rowIndex;
        delete [] store;
        rowIndex = new mkIndex[maxnnz];
        store = new Real[maxnnz];
    }

    if (ncols != B.Ncols()) {
        // re-allocate precisely required memory for column pointers
        ncols = B.Ncols();
        delete [] colPointer;
        colPointer = new mkIndex[ncols+1u];
    }

    // set number of rows
    nrows = B.Nrows();
    // copy all
    memory_xcopy(nnz, B.Store(), store);
    memcpy(rowIndex, B.RowIndex(), nnz*sizeof(mkIndex));
    memcpy(colPointer, B.ColPointer(), (ncols+1u)*sizeof(mkIndex));

    // return the result
    return *this;
}

// operator -
// return -A, where A is *this
// this operator invokes vector_scalar_operation() for copying elements
SparseMatrix SparseMatrix::operator-() const {
    // number of nonzero elements
    mkIndex nnz = colPointer[ncols]-colPointer[0];

    // allocate memory for row indices of output and copy the values
    mkIndex *rowIndex2 = new mkIndex[nnz];
    memcpy(rowIndex2, RowIndex(), nnz*sizeof(mkIndex));

    // allocate memory for column pointers of output and copy the values
    mkIndex *colPointer2 = new mkIndex[ncols+1u];
    memcpy(colPointer2, ColPointer(), (ncols+1u)*sizeof(mkIndex));

    // pointer for output
    Real *store2 = NULL;

    // elementwise multiply -1.0
    vector_scalar_operation(nnz, Store(), -1.0, store2, 2);

    // finally, return the result
    return SparseMatrix(nrows, ncols, store2, rowIndex2, colPointer2);
}


// operator *,/,*=,/= a real number
// return B with B(i,j)=A(i,j)*f for i,j=1,...,nrows, where A is *this
// this operator invokes vector_scalar_operation() to copy scaled nonzero elements
SparseMatrix SparseMatrix::operator*(Real f) const {
    // a special case: f == 0
    if (f == 0.0)
        return SparseMatrix(nrows, ncols);

    // number of nonzero elements
    mkIndex nnz = colPointer[ncols]-colPointer[0];

    // allocate memory for row indices of output and copy the values
    mkIndex *rowIndex2 = new mkIndex[nnz];
    memcpy(rowIndex2, RowIndex(), nnz*sizeof(mkIndex));

    // allocate memory for column pointers of output and copy the values
    mkIndex *colPointer2 = new mkIndex[ncols+1u];
    memcpy(colPointer2, ColPointer(), (ncols+1u)*sizeof(mkIndex));

    // pointer for output
    Real *store2 = NULL;
    // copy the scaled elements
    vector_scalar_operation(nnz, Store(), f, store2, 2);

    // return the result
    return SparseMatrix(nrows, ncols, store2, rowIndex2, colPointer2);
}

// scalar - sparse matrix multiplication
// return f*A
// this operator (eventually) invokes vector_scalar_operation() to copy scaled nonzero elements
SparseMatrix operator*(Real f, const SparseMatrix &A) {
    return A*f;
}

// return B with B(i,j)=A(i,j)/f for i,j=1,...,nrows, where A is *this
// this operator invokes vector_scalar_operation() to copy scaled nonzero elements
SparseMatrix SparseMatrix::operator/(Real f) const {
    // exception handling
    if (f == 0.0) {
        *Basics::err << "SparseMatrix::operator/=(Real): division-by-zero exception!" << endl;
        Basics::quit(1);
    }

    // number of nonzero elements
    mkIndex nnz = colPointer[ncols]-colPointer[0];

    // allocate memory for row indices of output and copy the values
    mkIndex *rowIndex2 = new mkIndex[nnz];
    memcpy(rowIndex2, RowIndex(), nnz*sizeof(mkIndex));

    // allocate memory for column pointers of output and copy the values
    mkIndex *colPointer2 = new mkIndex[ncols+1u];
    memcpy(colPointer2, ColPointer(), (ncols+1u)*sizeof(mkIndex));

    // pointer for output
    Real *store2 = NULL;
    // copy the scaled elements
    vector_scalar_operation(nnz, Store(), f, store2, 3);

    // return the result
    return SparseMatrix(nrows, ncols, store2, rowIndex2, colPointer2);
}

// return a const reference to A after computing A(i,j)*=f for i=1,...,nrows and j=1,...,ncols, where A is *this
// this operator invokes vector_scalar_operation() to scale nonzero elements
const SparseMatrix& SparseMatrix::operator*=(Real f) {
    // scale the elements
    mkIndex nnz = colPointer[ncols]-colPointer[0];
    vector_scalar_operation(nnz, Store(), f, store, 2);

    // return the result
    return *this;
}

// return a const reference to A after computing A(i,j)/=f for i=1,...,nrows and j=1,...,ncols, where A is *this
// this operator invokes vector_scalar_operation() to scale nonzero elements
const SparseMatrix& SparseMatrix::operator/=(Real f) {
    // exception handling
    if (f == 0.0) {
        *Basics::err << "SparseMatrix::operator/=(Real): division-by-zero exception!" << endl;
        Basics::quit(1);
    }

    // scale the elements
    mkIndex nnz = colPointer[ncols]-colPointer[0];
    vector_scalar_operation(nnz, Store(), f, store, 3);

    // return the result
    return *this;
}


// sparse matrix - matrix addition/subtraction
// return A+B, where A is *this
// this member function (eventually) invokes matrix_spmatrix_addition()
Matrix SparseMatrix::operator+(const Matrix &A) const {
    return A+(*this);
}

// return A-B, where A is *this
// this member function (eventually) invokes matrix_spmatrix_addition()
Matrix SparseMatrix::operator-(const Matrix &A) const {
    return -A+(*this);
}


// operator +,-,+=,-= a sparse matrix
// return A+B, where A is *this
// this member function invokes spmatrix_spmatrix_addition()
SparseMatrix SparseMatrix::operator+(const SparseMatrix &B) const {
    // exception handling
    if (nrows!=B.Nrows() || ncols!=B.Ncols()) {
        *Basics::err << "SparseMatrix::operator+(const SparseMatrix&): matrix dimensions must agree!" << endl;
        Basics::quit(1);
    }

    // pointers for output
    mkIndex *colPtr2 = NULL;
    mkIndex *rowIdx2 = NULL;
    Real *store2 = NULL;

    // compute C = A+B, where A is *this and C is a sparse matrix stored in CSC format in colPtr2[], rowIdx2[], store2[]
    spmatrix_spmatrix_addition(ncols, 1.0, Store(), RowIndex(), ColPointer(), 1.0, B.Store(), B.RowIndex(), B.ColPointer(), store2, rowIdx2, colPtr2);

    // return the result
    return SparseMatrix(nrows, ncols, store2, rowIdx2, colPtr2, Nnz()+B.Nnz());
}

// return A-B, where A is *this
// this member function invokes spmatrix_spmatrix_addition()
SparseMatrix SparseMatrix::operator-(const SparseMatrix &B) const {
    if (nrows!=B.Nrows() || ncols!=B.Ncols()) {
        *Basics::err << "SparseMatrix::operator-(const SparseMatrix&) const: matrix dimensions must agree!" << endl;
        Basics::quit(1);
    }

    // pointers for output
    mkIndex *colPtr2a = NULL;
    mkIndex *rowIdx2a = NULL;
    Real *store2a = NULL;

    // compute C = A-B, where A is *this and C is a sparse matrix stored in CSC format in colPtr2[], rowIdx2[], store2[]
    spmatrix_spmatrix_addition(ncols, 1.0, Store(), RowIndex(), ColPointer(), -1.0, B.Store(), B.RowIndex(), B.ColPointer(), store2a, rowIdx2a, colPtr2a);

    // return sparse matrix C = A - B
    return SparseMatrix(nrows, ncols, store2a, rowIdx2a, colPtr2a, Nnz()+B.Nnz());
}

// return a const reference to A after computing A = A+B, where A is *this
// this member function (eventually) invokes spmatrix_spmatrix_addition()
const SparseMatrix& SparseMatrix::operator+=(const SparseMatrix &B) {
    if (nrows!=B.Nrows() || ncols!=B.Ncols()) {
        *Basics::err << "SparseMatrix::operator+=(const SparseMatrix&): matrix dimensions must agree!" << endl;
        Basics::quit(1);
    }

    // compute A = A+B, where A is *this
    *this = *this + B;

    // return the result
    return *this;
}

// return a const reference to A after computing A = A-B, where A is *this
// this member function (eventually) invokes spmatrix_spmatrix_addition()
const SparseMatrix& SparseMatrix::operator-=(const SparseMatrix &B) {
    if (nrows!=B.Nrows() || ncols!=B.Ncols()) {
        *Basics::err << "SparseMatrix::operator-=(const SparseMatrix&): matrix dimensions must agree!" << endl;
        Basics::quit(1);
    }

    // compute A = A-B, where A is *this
    *this = *this - B;

    // return the result
    return *this;
}


// sparse matrix - vector multiplication
// return A*v, where A is *this
// this operator invokes spmatrix_vector_multiplication()
Vector SparseMatrix::operator*(const Vector &v) const {
    if (ncols != v.Length()) {
        *Basics::err << "SparseMatrix::operator*(const Vector&) const: inner matrix-vector dimensions must agree!" << endl;
        Basics::quit(1);
    }

    // pointer for output
    Real *sw = NULL;

    // compute w = A*v, where A is *this
    spmatrix_vector_multiplication(nrows, ncols, store, rowIndex, colPointer, v.Store(), sw);

    // return the result
    return Vector(nrows, sw);
}

// sparse matrix - matrix multiplication
// return A*B, where A is *this
// this operator invokes spmatrix_matrix_multiplication()
Matrix SparseMatrix::operator*(const Matrix &B) const {
    if (ncols != B.Nrows()) {
        *Basics::err << "SparseMatrix::operator*(const Matrix&) const: inner matrix dimensions must agree!" << endl;
        Basics::quit(1);
    }

    // pointer for output
    Real *sc = NULL;

    // compute C = A * B, where A is *this
    spmatrix_matrix_multiplication(nrows, ncols, B.Ncols(), Store(), RowIndex(), ColPointer(), B.Store(), sc);

    // return the result
    return Matrix(nrows, B.Ncols(), sc);
}


// the member function permuteColumns permutes columns, such that
//    B(i,perm[j-1])   = A(i,j) for i=1,...,nr and j=1,...,nc (indexFrom1==true), or
//    B(i,perm[j-1]+1) = A(i,j) for i=1,...,nr and j=1,...,nc (indexFrom1==false), where
//    A is *this and B is the returned sparse matrix
// this member function invokes permute_csc_columns()
// this member function resembles B(:,perm) = A in OCTAVE/MATLAB
SparseMatrix SparseMatrix::permuteColumns(const mkIndex *perm, bool indexFrom1) const {
    // pointers for output
    Real *store2 = NULL;
    mkIndex *rowIdx2 = NULL;
    mkIndex *colPtr2 = NULL;

    // compute the permutation
    permute_csc_columns(ncols, perm, store, rowIndex, colPointer, store2, rowIdx2, colPtr2, indexFrom1);

    // return the result
    return SparseMatrix(nrows, ncols, store2, rowIdx2, colPtr2);
}

// the member function permuteRows permutes rows, such that
//    B(perm[i-1],j)   = A(i,j) for i=1,...,nr and j=1,...,nc (indexFrom1==true), or
//    B(perm[i-1]+1,j) = A(i,j) for i=1,...,nr and j=1,...,nc (indexFrom1==false), where
//    A is *this and B is the return sparse matrix
// the indices perm[0,...,nrows-1] can be 0,1,...,nrows-1 (indexFrom1==false) or 1,2,...,nrows (indexFrom1==true)
// this member function invokes permute_csc_rows() and then sort_csc_elements_in_each_column(),
//    so that the return sparse matrix has elements in each column sorted with respect to row indices
// this member function resembles B(perm,:) = A in OCTAVE/MATLAB
SparseMatrix SparseMatrix::permuteRows(const mkIndex *perm, bool indexFrom1) const {
    // pointers for output
    Real *store2 = NULL;
    mkIndex *rowIdx2 = NULL;
    mkIndex *colPtr2 = NULL;

    // compute the permutation
    permute_csc_rows(ncols, perm, store, rowIndex, colPointer, store2, rowIdx2, colPtr2, indexFrom1);

    // after permutation of rows, the elements in each column may not be sorted with respect to row indices
    // so we sort the elements in each column
    sort_csc_elements_in_each_column(nrows, ncols, store2, rowIdx2, colPtr2);
    // now the elements in each column is sorted with respect to row indices

    // return the result
    return SparseMatrix(nrows, ncols, store2, rowIdx2, colPtr2);
}

// the member functionpermuteRowsAndColumns permutes rows and columns, such that
//    B(rperm[i-1],cperm[j-1])     = A(i,j) for i=1,...,nr and j=1,...,nc (if indexFrom1==true), or
//    B(rperm[i-1]+1,cperm[j-1]+1) = A(i,j) for i=1,...,nr and j=1,...,nc (if indexFrom1==false), where
//    A is the input matrix *this and B is the returned output matrix
// if cperm==NULL, then symmetric permutation will be applied (i.e., *rperm will also be used for *cperm)
// the indices rperm[0,...,nrows-1] range over 0,1,...,nrows-1 (if indexFrom1==false) or over 1,2,...,nrows (if indexFrom1==true)
// the indices cperm[0,...,ncols-1] range over 0,1,...,ncols-1 (if indexFrom1==false) or over 1,2,...,ncols (if indexFrom1==true)
// this member function resembles B(rperm,cperm) = A in OCTAVE/MATLAB
SparseMatrix SparseMatrix::permuteRowsAndColumns(const mkIndex *rperm, const mkIndex *cperm, bool indexFrom1) const {
    if (cperm == NULL) {
        if (nrows != ncols) {
            *Basics::err << "SparseMatrix::permuteRowsAndColumns(const mkIndex*, const mkIndex*, bool) const: symmetric permutation requires the sparse matrix being square!" << endl;
            Basics::quit(1);
        }
        return permuteRows(rperm, indexFrom1).permuteColumns(rperm, indexFrom1);
    }
    return permuteRows(rperm, indexFrom1).permuteColumns(cperm, indexFrom1);
}


#ifdef USE_MEX
// convert mxArray to SparseMatrix
SparseMatrix SparseMatrix::mxArray2SparseMatrix(const mxArray *mA) {
    if (!mxIsNumeric(mA))
        mexErrMsgTxt("mxArray2SpasreMatrix(const mxArray*): the input mxArray* contains non-numeric data!");
    if (!mxIsSparse(mA))
        mexErrMsgTxt("mxArray2SparseMatrix(const mxArray*): the input mxArray* is not a sparse matrix!");

    // OCTAVE/MATLAB: [nrows, ncols] = size(A);
    mkIndex nrows = (mkIndex)mxGetM(mA);
    mkIndex ncols = (mkIndex)mxGetN(mA);

    // OCTAVE/MATLAB: [I, J, V] = find(A);
    mwIndex *ii = mxGetIr(mA);
    mwIndex *jj = mxGetJc(mA);

    // allocate memory for indices for the output
    mkIndex nnz = jj[ncols] - jj[0];
    mkIndex *rowIndex = NULL;
    if (nnz)
        rowIndex = new mkIndex[nnz];
    mkIndex *colPointer = new mkIndex[ncols+1u];

    // copy the indices
    mkIndex z = nnz;
    mkIndex *ri = rowIndex;
    mwIndex *ij = ii;
    while (z--)
        *ri++ = (mkIndex)(*ij++);
    mkIndex c = ncols+1u;
    mkIndex *cp = colPointer;
    ij = jj;
    while (c--)
        *cp++ = (mkIndex)(*ij++);

    // copy the numerical data
    Real *store = NULL;
    if (nnz)
        store = new Real[nnz];
    if (mxIsDouble(mA)) {
        double *vv = (double *)mxGetPr(mA);
        #ifdef USE_SINGLE
            mexWarnMsgTxt("mxArray2SparseMatrix(const mxArray*): the input mxArray* contains double precision floating-point data, but the output will be in single precision!");
            z = nnz;
            Real *s = store;
            while (z--)
                *s++ = (Real)(*vv++);
        #else
            memory_xcopy(nnz, vv, store);
        #endif
    }
    else if (mxIsSingle(mA)) {
        float *vv = (float *)mxGetData(mA);
        #ifndef USE_SINGLE
            mexWarnMsgTxt("mxArray2SparseMatrix(const mxArray*): the input mxArray* contains single precision floating-point data, but the output will be in double precision!");
            z = nnz;
            Real *s = store;
            while (z--)
                *s++ = (Real)(*vv++);
        #else
            memory_xcopy(nnz, vv, store);
        #endif
    }
    else {
        mexErrMsgTxt("mxArray2SparseMatrix(const mxArray*): the input mxArray* contains integer-type data; cannot convert it!");
    }
    if (mxIsComplex(mA))
        mexWarnMsgTxt("mxArray2SparseMatrix(const mxArray*): the input mxArray* is a complex sparse matrix; the imaginary part is ignored!");

    // return the result
    return SparseMatrix(nrows, ncols, store, rowIndex, colPointer);
}

// convert SparseMatrix to mxArray
mxArray *SparseMatrix::SparseMatrix2mxArray(const SparseMatrix &A) {
    #ifdef USE_SINGLE
        mexWarnMsgTxt("SparseMatrix2MxArray(const SparseMatrix &): the input SparseMatrix is in single precision; the converted mxArray* will be in double precision!");
    #endif

    // create the mxArray and get its pointers
    mxArray *mA = mxCreateSparse((mwSize)A.Nrows(), (mwSize)A.Ncols(), (mwSize)A.Nnz(), mxREAL);
    double  *sa = mxGetPr(mA);
    mwIndex *ii = mxGetIr(mA);
    mwIndex *jj = mxGetJc(mA);

    // copy the data
    const mkIndex *ri = A.RowIndex();
    const mkIndex *cp = A.ColPointer();
    const Real *str = A.Store();
    mkIndex z = A.Nnz();
    while (z--) {
        *ii++ = (mwIndex)(*ri++);
        *sa++ = (Real)(*str++);
    }
    mkIndex c = A.Ncols()+1;
    while (c--)
        *jj++ = (mwIndex)(*cp++);

    // return the result
    return mA;
}
#endif  // of #ifdef USE_MEX


#ifdef USE_NAMESPACE
}  // end of namespace NEWMAT
#endif
