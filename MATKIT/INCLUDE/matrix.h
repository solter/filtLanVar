#ifndef MATRIX_H
#define MATRIX_H
/*! \file matrix.h
 *  The class Matrix defines a general matrix represented by
 *  an array Matrix::store[] (with elements of the matrix stored columnwise) and
 *  the size Matrix::nrows by Matrix::ncols.
 *
 *  To be precise, The (<em>i,j</em>) entry of the matrix is stored as
 *  Matrix::store[(<em>i-1</em>)+(<em>j-1</em>)*<em>nrows</em>]
 *  for <em>i=1,...,nrows</em> and <em>j=1,...,ncols</em>.
 */

#include <iostream>  // for ostream, etc.
#include "matkitdef.h"


#ifdef USE_NAMESPACE
namespace MATKIT {
#endif


class Vector;
class SymmetricMatrix;
class SparseMatrix;

/*! \brief Class Matrix defines a generic matrix and the arithmetic operations on it,
 *         and also providing some static functions for matrix generation.
 */
class Matrix {
    // friend classes and functions
    #ifdef USE_MWBLAS
    // MWBLAS routines dgemv (double precision) and sgemv (single precision) for vector - matrix multiplication demand non-const access to the matrix but indeed only const access is required
    friend Vector Vector::operator*(const Matrix &A) const;
    friend const Vector &Vector::operator*=(const Matrix &A);
    #endif

protected:
    //! The matrix <em>*this</em> is stored columnwise in <em>store</em>[].
    /*! To be precise, let <em>A</em> be the matrix <em>*this</em>. Then
     *  <em>A</em>(<em>i,j</em>) is stored in <em>store</em>[(<em>i-1</em>)+(<em>j-1</em>)<em>*nrows</em>]
     *  for <em>i=1,...,nrows</em> and <em>j=1,...,ncols</em>.
     */
    // the matrix *this is stored columnwise in store[]
    // to be precise, let A be the matrix *this
    // then A(i,j) is stored in store[(i-1)+(j-1)*nrows] for i=1,...,nrows and j=1,...,ncols
    Real *store;
    // number of rows and number of columns of *this
    //! Number of rows of <em>*this</em>.
    mkIndex nrows;
    //! Number of columns of <em>*this</em>.
    mkIndex ncols;
    //! Maximum number of columns allowed in the allocated memory.
    /*! While memory of size <em>nrows*ncols</em> for storing the matrix is required,
     *  more memory (with the address pointed to by <em>store</em>) can be allocated.
     *  The current memory allows <em>nrows*maxncols</em> <b>Real</b> data (if <em>maxncols>=ncols</em>).
     *
     *  \remark  The reason to allocate more memory is to reduce the frequency of memory reallocation
     *           by Matrix::appendColumn().
     */
    // maximum number of columns allowed in the allocated memory
    mkIndex maxncols;

public:
    //! When the allocated memory is insufficient, this is the factor by which the memory is expanded (default <em>1.2</em>).
    /*! \remark  This is a parameter used by Matrix::appendColumn().
     *
     *  \see  Matrix::setMemoryExpansionFactor().
     */
    // when the allocated memory is insufficient, this is the factor by which the memory is expanded (default 1.2)
    static Real memoryExpansionFactor;

    // some static functions
    //! Return an <em>n</em>-by-<em>n</em> matrix of zeros.
    /*! \remark  This static function resembles <b>zeros</b>(<em>n</em>) in <b>OCTAVE</b>/<b>MATLAB</b>.
     */
    // return an n-by-n matrix of zeros
    // this static function resembles zeros(n) in OCTAVE/MATLAB
    static Matrix zeros(mkIndex n) { return zeros(n,n); }
    //! Return an <em>nr</em>-by-<em>nc</em> matrix of zeros.
    /*! \remark  This static function resembles <b>zeros</b>(<em>nr</em>,<em>nc</em>) in <b>OCTAVE</b>/<b>MATLAB</b>.
     */
    // return an nr-by-nc matrix of zeros
    // this static function resembles zeros(nr,nc) in OCTAVE/MATLAB
    static Matrix zeros(mkIndex nr, mkIndex nc);
    //! Return an <em>n</em>-by-<em>n</em> matrix of ones.
    /*! \remark  This static function resembles <b>ones</b>(<em>n</em>) in <b>OCTAVE</b>/<b>MATLAB</b>.
     */
    // return an n-by-n matrix of ones
    // this static function resembles ones(n) in OCTAVE/MATLAB
    static Matrix ones(mkIndex n) { return ones(n,n); }
    //! Return an <em>nr</em>-by-<em>nc</em> matrix of ones.
    /*! \remark  This static function resembles <b>ones</b>(<em>nr</em>,<em>nc</em>) in <b>OCTAVE</b>/<b>MATLAB</b>.
     */
    // return an nr-by-nc matrix of ones
    // this static function resembles ones(nr,nc) in OCTAVE/MATLAB
    static Matrix ones(mkIndex nr, mkIndex nc);
    //! Return an <em>n</em>-by-<em>n</em> matrix with ones on the main diagonal and zeros elsewhere.
    /*! \remark  This static function resembles <b>eye</b>(<em>n</em>) in <b>OCTAVE</b>/<b>MATLAB</b>.
     */
    // return an n-by-n matrix with ones on the main diagonal and zeros elsewhere
    // this static function resembles eye(n) in OCTAVE/MATLAB
    static Matrix eye(mkIndex n) { return eye(n,n); }
    //! Return an <em>nr</em>-by-<em>nc</em> matrix with ones on the main diagonal and zeros elsewhere.
    /*! \remark  This static function resembles <b>eye</b>(<em>nr</em>,<em>nc</em>) in <b>OCTAVE</b>/<b>MATLAB</b>.
     */
    // return an nr-by-nc matrix with ones on the main diagonal and zeros elsewhere
    // this static function resembles eye(nr,nc) in OCTAVE/MATLAB
    static Matrix eye(mkIndex nr, mkIndex nc);
    //! Return an <em>n</em>-by-<em>n</em> matrix containing pseudo-random values drawn from a uniform distribution in the interval [<em>a,b</em>].
    /*! \remark  This static function resembles <b>rand</b>(<em>n</em>) in <b>OCTAVE</b>/<b>MATLAB</b>.
     */
    // return an n-by-n matrix containing pseudo-random values drawn from a uniform distribution in the interval [a,b]
    // this static function resembles rand(n) in OCTAVE/MATLAB
    static Matrix random(mkIndex n, Real a=0.0, Real b=1.0) { return random(n,n,a,b); }
    //! Return an <em>nr</em>-by-<em>nc</em> matrix containing pseudo-random values drawn from a uniform distribution on the interval <em>[a,b]</em>.
    /*! \remark  This static function resembles <b>rand</b>(<em>nr</em>,<em>nc</em>) in <b>OCTAVE</b>/<b>MATLAB</b>.
     */
    // return an nr-by-nc matrix containing pseudo-random values drawn from a uniform distribution on the interval [a,b]
    // this static function resembles rand(nr,nc) in OCTAVE/MATLAB
    static Matrix random(mkIndex nr, mkIndex nc, Real a=0.0, Real b=1.0);


    //! Horizontal concatenation of two matrices <em>A</em> and <em>B</em>.
    /*! <em>A</em> and <em>B</em> must have the same number of rows.
     *  \remark  This static function invokes horizontal_concatenate().
     *  \remark  This static function resembles [<em>A</em>,<em>B</em>] or equivalently <b>horzcat</b>(<em>A</em>,<em>B</em>) in <b>OCTAVE</b>/<b>MATLAB</b>.
     */
    // horizontal concatenation of two matrices A and B
    // A and B must have the same number of rows
    // this static function invokes horizontal_concatenate()
    // this static function resembles [A,B] or equivalently horzcat(A,B) in OCTAVE/MATLAB
    static Matrix horizontalConcatenate(const Matrix &A, const Matrix &B);

    //! Vertical concatenation of two matrices <em>A</em> and <em>B</em>.
    /*! <em>A</em> and <em>B</em> must have the same number of columns.
     *  \remark  This static function invokes vertical_concatenate().
     *  \remark  This static function resembles [<em>A</em>;<em>B</em>] or equivalently <b>vertcat</b>(<em>A</em>,<em>B</em>) in <b>OCTAVE</b>/<b>MATLAB</b>.
     */
    // vertical concatenation of two matrices A and B
    // A and B must have the same number of columns
    // this static function invokes vertical_concatenate()
    // this static function resembles [A;B] or equivalently vertcat(A,B) in OCTAVE/MATLAB
    static Matrix verticalConcatenate(const Matrix &A, const Matrix &B);

    //! Sum of two matrices A and B, padding zeros if they have different dimensions.
    /*! \remark  This static function invokes matrix_xsum().
     */
    // sum of two matrices A and B, padding zeros if they have different dimensions
    // this static function invokes matrix_xsum()
    static Matrix xsum(const Matrix &A, const Matrix &B);


    // constructors
    //! A constructor for an empty matrix.
    // a constructor for an empty matrix
    Matrix();

    //! A constructor for an <em>nr</em>-by-<em>nc</em> matrix of zeros.
    /*! <ul>
     *  <li>  If <em>maxnc<=nc</em>, then memory of size <em>nr*nc</em> will be allocated.
     *  <li>  If <em>maxnc>nc</em>, then memory of size <em>nr*maxnc</em> will be allocated.
     *  </ul>
     */
    // a constructor for an nr-by-nc matrix of zeros
    // if maxnc<=nc, then memory of size nr*nc    will be allocated
    // if maxnc>nc,  then memory of size nr*maxnc will be allocated
    Matrix(mkIndex nr, mkIndex nc, mkIndex maxnc=0);
    //! A constructor for a matrix <em>A</em> of size <em>nr</em>-by-<em>nc</em> with elements stored columnwise in <em>store0</em>[].
    /*! <ul>
     *  <li>  If <em>maxnc<=nc</em>, then it is assumed that memory of size at least <em>nr*nc</em> has been allocated, with the address pointed to by <em>store0</em>.
     *  <li>  If <em>maxnc>nc</em>, then it is assumed that memory of size at least <em>nr*maxnc</em> has been allocated, with the address pointed to by <em>store0</em>.
     *  </ul>
     *  \remark It performs a shallow copy of <em>store0</em>[] and the memory will be freed by the destructor (i.e. `<em>delete</em> [] <em>store0</em>').
     */
    // a constructor for an nr-by-nc matrix of with elements stored columnwise in store0[]
    // if maxnc<=nc, then it is assumed that memory of size at least nr*nc    has been allocated, with the address pointed to by store0
    // if maxnc>nc,  then it is assumed that memory of size at least nr*maxnc has been allocated, with the address pointed to by store0
    // note that it performs a shallow copy of store0 and the memory will be freed by the destructor (i.e. `delete [] store0')
    Matrix(mkIndex nr, mkIndex nc, Real *store0, mkIndex maxnc=0);

    //! A copy constructor.
    /*! \remark  This constructor invokes memory_xcopy().
     */
    // a copy constructor
    // this constructor invokes memory_xcopy()
    Matrix(const Matrix &B);

    //! A destructor. Memory allocated, with the address pointed to by Matrix::store, will be freed.
    // a destructor; memory allocated, with the address pointed to by Matrix::store, will be freed
    virtual ~Matrix();


    // get the information of *this
    //! The number of rows of <em>*this</em>.
    // the number of rows of *this
    mkIndex Nrows() const { return nrows; }
    //! The number of columns of <em>*this</em>.
    // the number of columns of *this
    mkIndex Ncols() const { return ncols; }
    //! The number of maximum number of columns allowed in the allocated memory.
    // the number of maximum number of columns allowed in the allocated memory
    mkIndex Maxncols() const { return maxncols; }
    //! A <b>const</b> pointer to the storage for elements of <em>*this</em>.
    // a const pointer to the storage for elements of *this
    const Real *Store() const { return store; }


    //! Convert <em>*this</em> to a sparse matrix.
    /*! \remark  This member function invokes full_to_sparse().
     *  \remark  This member function resembles <b>sparse</b>(<em>A</em>) in <b>OCTAVE</b>/<b>MATLAB</b>, where <em>A</em> is <em>*this</em>.
     */
    // convert *this to a sparse matrix
    // this member function invokes full_to_sparse()
    // this member function resembles sparse(A) in OCTAVE/MATLAB, where A is *this
    SparseMatrix toSparse() const;

    //! Resize <em>*this</em> to be of size <em>nr</em>-by-<em>nc</em>.
    /*! <ul>
     *  <li>  If <em>lazyMalloc</em>==<b>false</b>, then memory will be allocated whenever the size of allocated memory is different from required.
     *  <li>  If <em>lazyMalloc</em>==<b>true</b>, then memory will be reallocated only when the allocated memory is insufficient.
     *  </ul>
     */
    // resize *this to be of size nr-by-nc
    // if lazyMalloc == false, then memory will be allocated whenever the size of allocated memory is different from required
    // if lazyMalloc == true, then memory will be allocated only when the allocated memory is insufficient
    void resize(mkIndex nr, mkIndex nc, bool lazyMalloc=true);

    //! Set <em>*this</em> be an <em>nr</em>-by-<em>nc</em> matrix stored in <em>store0</em>[].
    /*  <ul>
     *  <li>  If <em>maxnc<=nc</em>, then it is assumed that memory of size at least <em>nr*nc</em> has been allocated, with the address pointed to by <em>store0</em>.
     *  <li>  If <em>maxnc>nc</em>, then it is assumed that memory of size at least <em>nr*maxnc</em> has been allocated, with the address pointed to by <em>store0</em>.
     *  </ul>
     *  \remark It performs a shallow copy of <em>store0</em>[] and the memory will be freed by the destructor (i.e. `<em>delete</em> [] <em>store0</em>').
     */
    // set *this be an nr-by-nc matrix stored in store0[]
    // if maxnc<=nc, then it is assumed that memory of size at least nr*nc    has been allocated, with the address pointed to by store0
    // if maxnc>nc,  then it is assumed that memory of size at least nr*maxnc has been allocated, with the address pointed to by store0
    // it performs a shallow copy of store0[] and the memory will be freed by the destructor (i.e. `delete [] store0')
    void set(mkIndex nr, mkIndex nc, Real *store0, mkIndex maxnc=0);


    //! Append a column formed by <em>vs</em>[<em>0,...,nrows-1</em>] to <em>*this</em>.
    /*! \remark  This member function invokes memory_xcopy().
     *  \remark  Allocated memory, if insufficient, will be expanded by a factor of Matrix::memoryExpansionFactor.
     */
    // append a column formed by vs[0,...,nrows-1] to *this
    // this member function invokes memory_xcopy()
    // allocated memory, if insufficient, will be expanded by a factor of Matrix::memoryExpansionFactor
    void appendColumn(const Real *vs);
    //! Append a column formed by <em>v</em> to <em>*this</em>.
    /*! \remark  This member function (eventually) invokes memory_xcopy().
     *  \remark  Allocated memory, if insufficient, will be expanded by a factor of Matrix::memoryExpansionFactor.
     */
    // append a column formed by v to *this
    // this member function (eventually) invokes memory_xcopy()
    // allocated memory, if insufficient, will be expanded by a factor of Matrix::memoryExpansionFactor
    void appendColumn(const Vector &v);

    //! Set the memory expansion factor.
    /*! The memory expansion factor determines how much memory should be reallocated
     *  when the allocated memory is insufficient to store the matrix in Matrix::appendColumn().
     *
     *  \see  Matrix::memoryExpansionFactor.
     */
    // set memory expansion factor
    // the memory expansion factor determines how much memory should be reallocated
    // when the allocated memory is insufficient to store the matrix in Matrix::appendColumn()
    // the default value of memoryExpansionFactor is 1.2, shown in matrix.cpp
    // see also Matrix::memoryExpansionFactor
    static void setMemoryExpansionFactor(Real newExpansionFactor);


    //! A vector formed by column <em>j</em> of <em>*this</em>.
    /*! \remark  This member function invokes memory_xcopy().
     *  \remark  This member function resembles <em>A</em>(:,<em>j</em>) in <b>OCTAVE</b>/<b>MATLAB</b>, where <em>A</em> is <em>*this</em>.
     *  \see  Matrix::row(), Matrix::rows(), Matrix::columns(), Matrix::subMatrix().
     */
    // a vector formed by column j of *this
    // this member function invokes memory_xcopy()
    // this member function resembles A(:,j) in OCTAVE/MATLAB, where A is *this
    Vector column(mkIndex j) const;
    //! A vector formed by row <em>j</em> of <em>*this</em>.
    /*! \remark  This member function invokes memory_xcopy().
     *  \remark  This member function resembles <em>A</em>(<em>i</em>,:) in <b>OCTAVE</b>/<b>MATLAB</b>, where <em>A</em> is <em>*this</em>.
     *  \see  Matrix::column(), Matrix::rows(), Matrix::columns(), Matrix::subMatrix().
     */
    // a vector formed by row i of *this
    // this member function invokes memory_xcopy()
    // this member function resembles A(i,:) in OCTAVE/MATLAB, where A is *this
    Vector row(mkIndex i) const;

    //! A submatrix formed by columns <em>j1,...,j2</em> of <em>*this</em>.
    /*! \remark  This member function invokes submatrix_of_matrix().
     *  \remark  This member function resembles <em>A</em>(:,<em>j1:j2</em>) in <b>OCTAVE</b>/<b>MATLAB</b>, where <em>A</em> is <em>*this</em>.
     *  \see  Matrix::row(), Matrix::column(), Matrix::rows(), Matrix::subMatrix().
     */
    // a submatrix formed by columns j1,...,j2 of *this
    // this member function invokes submatrix_of_matrix()
    // this member function resembles A(:,j1:j2) in OCTAVE/MATLAB, where A is *this
    Matrix columns(mkIndex j1, mkIndex j2) const;
    //! A submatrix formed by rows <em>i1,...,i2</em> of <em>*this</em>.
    /*! \remark  This member function invokes submatrix_of_matrix().
     *  \remark  This member function resembles <em>A</em>(<em>i1:i2</em>,:) in <b>OCTAVE</b>/<b>MATLAB</b>, where <em>A</em> is <em>*this</em>.
     *  \see  Matrix::row(), Matrix::column(), Matrix::columns(), Matrix::subMatrix().
     */
    // a submatrix formed by rows i1,...,i2 of *this
    // this member function invokes submatrix_of_matrix()
    // this member function resembles A(i1:i2,:) in OCTAVE/MATLAB, where A is *this
    Matrix rows(mkIndex i1, mkIndex i2) const;
    //! A submatrix formed by elements in rows <em>i1,...,i2</em> and columns <em>j1,...,j2</em> of <em>*this</em>.
    /*! \remark  This member function invokes submatrix_of_matrix().
     *  \remark  This member function resembles <em>A</em>(<em>i1:i2,j1:j2</em>) in <b>OCTAVE</b>/<b>MATLAB</b>, where <em>A</em> is <em>*this</em>
     *  \see  Matrix::row(), Matrix::column(), Matrix::rows(), Matrix::columns().
     */
    // a submatrix formed by elements in rows i1,...,i2 and columns j1,...,j2 of *this
    // this member function invokes submatrix_of_matrix()
    // this member function resembles A(i1:i2,j1:j2) in OCTAVE/MATLAB, where A is *this
    Matrix subMatrix(mkIndex i1, mkIndex i2, mkIndex j1, mkIndex j2) const;


    //! A submatrix formed by columns <em>c</em>[<em>0</em>],<em>c</em>[<em>1</em>],...,<em>c</em>[<em>len-1</em>] of <em>*this</em>.
    /*! \remark  This member function resembles <em>A</em>(:,<em>c</em>) in <b>OCTAVE</b>/<b>MATLAB</b>, where <em>A</em> is <em>*this</em>.
     *  \see  Matrix::row(), Matrix::column(), Matrix::rows(), Matrix::subMatrix().
     */
    // a submatrix formed by columns c[0],c[1],...,c[len-1] of *this
    // this member function invokes memory_xcopy() iteratively
    // this member function resembles A(:,c) in OCTAVE/MATLAB, where A is *this
    Matrix columns(mkIndex len, const mkIndex c[], bool indexFrom1=true) const;
    //! A submatrix formed by rows <em>r</em>[<em>0</em>],<em>r</em>[<em>1</em>],...,<em>r</em>[<em>len-1</em>] of <em>*this</em>.
    /*! \remark This member function resembles <em>A</em>(<em>r</em>,:) in <b>OCTAVE</b>/<b>MATLAB</b>, where <em>A</em> is <em>*this</em>.
     *  \see    Matrix::row(), Matrix::column(), Matrix::columns(), Matrix::subMatrix().
     */
    // a submatrix formed by rows r[0],r[1],...,r[len-1] of *this
    // this member function resembles A(r,:) in OCTAVE/MATLAB, where A is *this
    Matrix rows(mkIndex len, const mkIndex r[], bool indexFrom1=true) const;
    //! A submatrix formed by elements in rows <em>r</em>[<em>0</em>]<em>,r</em>[<em>1</em>],...,<em>r</em>[<em>rlen-1</em>] and columns <em>c</em>[<em>0</em>]<em>,c</em>[<em>1</em>],...,<em>c</em>[<em>clen-1</em>] of <em>*this</em>.
    /*! \remark  This member function resembles <em>A</em>(<em>r,c</em>) in <b>OCTAVE</b>/<b>MATLAB</b>, where <em>A</em> is <em>*this</em>.
     *  \see  Matrix::row(), Matrix::column(), Matrix::rows(), Matrix::columns().
     */
    // a submatrix formed by elements in rows r[0],r[1],...,r[rlen-1] and columns c[0],c[1],...,c[clen-1] of *this
    // this member function resembles A(r,c) in OCTAVE/MATLAB, where A is *this
    Matrix subMatrix(mkIndex rlen, const mkIndex r[], mkIndex clen, const mkIndex c[], bool indexFrom1=true) const;


    //! A vector formed by the elements in the <em>k</em>th diagonal of <em>*this</em> in order.
    /*! \remark  This member function invokes matrix_diag().
     *  \remark  This member function resembles <b>diag</b>(<em>A</em>,<em>k</em>) in <b>OCTAVE</b>/<b>MATLAB</b>.
     */
    // a vector formed by the elements in the k-th diagonal of *this in order
    // this member function invokes matrix_diag()
    // this member function resembles diag(A,k) in OCTAVE/MATLAB
    Vector diag(mkSignedIndex k=0) const;
    //! The transpose of <em>*this</em>.
    /*! \remark  This member function invokes matrix_transpose().
     *  \remark  This member function resembles <em>A</em>' or equivalently <b>transpose</b>(<em>A</em>) in <b>OCTAVE</b>/<b>MATLAB</b>,
     *           where <em>A</em> is <em>*this</em>.
     */
    // the transpose of *this
    // this member function invokes matrix_transpose()
    // this member function resembles A' or equivalently transpose(A) in OCTAVE/MATLAB, where A is *this
    Matrix transpose() const;


    // norms
    //! Frobenius norm of <em>*this</em>.
    // Frobenius norm of *this
    Real normFrobenius() const;
    //! 1-norm of <em>*this</em>.
    // 1-norm of *this
    Real norm1() const;
    //! \f$\infty\f$-norm of <em>*this</em>.
    // infinity-norm of *this
    Real normInfinity() const;
    //! Sum of squared elements of <em>*this</em>.
    // sum of squared elements of *this
    Real sumSquare() const;


    // access A(i,j)
    //! Return a reference to <em>A</em>(<em>i,j</em>), i.e. can be used to write <em>A</em>(<em>i,j</em>), where A is <em>*this</em>.
    // return a reference to A(i,j), i.e. can be used to write A(i,j), where A is *this
    virtual Real& operator()(mkIndex i, mkIndex j);
    //! Read <em>A</em>(<em>i,j</em>) via a <b>const</b> reference or pointer to <em>A</em>, where <em>A</em> is <em>*this</em>.
    // read A(i,j) via a const reference or pointer to A, where A is *this
    Real operator()(mkIndex i, mkIndex j) const;


    // operators ==, !=
    //! Return <b>true</b> if <em>A==B</em> mathematically, and otherwise return <b>false</b>, where <em>A</em> is <em>*this</em>.
    /*! \remark  This operator invokes are_elements_equal_elements().
     */
    // return true if A == B mathematically, and otherwise return false, where A is *this
    // this operator invokes are_elements_equal_elements()
    bool operator==(const Matrix &B) const {
        if (nrows != B.nrows || ncols != B.ncols)
            return false;
        return are_elements_equal_elements(nrows*ncols, Store(), B.Store());
    }
    //! Return <b>false</b> if <em>A==B</em> mathematically, and otherwise return <b>true</b>, where <em>A</em> is <em>*this</em>.
    /*! \remark  This operator (eventually) invokes are_elements_equal_elements().
     */
    // return false if A == B mathematically, and otherwise return true, where A is *this
    // this operator (eventually) invokes are_elements_equal_elements()
    bool operator!=(const Matrix &B) const {
        return !(*this==B);
    }


    // operator =
    //! Return a <b>const</b> reference to <em>A</em> after setting <em>A</em>(<em>i,j</em>)=<em>f</em> for <em>i=1,...,nrows</em> and <em>j=1,...,ncols</em>, where <em>A</em> is <em>*this</em>.
    // return a const reference to A after setting A(i,j)=f for i=1,...,nrows and j=1,...,ncols, where A is *this
    const Matrix& operator=(Real f);
    //! Return a <b>const</b> reference to <em>A</em> after setting <em>A=B</em>, where <em>A</em> is <em>*this</em>.
    /*! This operator invokes memory_xcopy().
     */
    // return a const reference to A after setting A=B, where A is *this
    // this operator invokes memory_xcopy()
    const Matrix& operator=(const Matrix &B);


    // operator -
    //! Return -<em>A</em>, where <em>A</em> is <em>*this</em>.
    /*! \remark  This operator invokes vector_scalar_operation().
     */
    // return -A, where A is *this
    // this operator invokes vector_scalar_operation()
    Matrix operator-() const;


    // operators +,-,*,/,+=,-=,*=,/= a real number
    //! Return <em>B</em> with <em>B</em>(<em>i,j</em>)=<em>A</em>(<em>i,j</em>)+<em>f</em> for <em>i=1,...,nrows</em> and <em>j=1,...,ncols</em>, where <em>A</em> is <em>*this</em>.
    /*! \remark  This operator invokes vector_scalar_operation().
     */
    // return B with B(i,j)=A(i,j)+f for i=1,...,nrows and j=1,...,ncols, where A is *this
    // this operator invokes vector_scalar_operation()
    Matrix operator+(Real f) const;
    //! Return <em>B</em> with <em>B</em>(<em>i,j</em>)=<em>A</em>(<em>i,j</em>)-<em>f</em> for <em>i=1,...,nrows</em> and <em>j=1,...,ncols</em>, where <em>A</em> is <em>*this</em>.
    /*! \remark  This operator invokes vector_scalar_operation().
     */
    // return B with B(i,j)=A(i,j)-f for i=1,...,nrows and j=1,...,ncols, where A is *this
    // this operator invokes vector_scalar_operation()
    Matrix operator-(Real f) const;
    //! Return <em>B</em> with <em>B</em>(<em>i,j</em>)=<em>A</em>(<em>i,j</em>)*<em>f</em> for <em>i=1,...,nrows</em> and <em>j=1,...,ncols</em>, where <em>A</em> is <em>*this</em>.
    /*! \remark  This operator invokes vector_scalar_operation().
     */
    // return B with B(i,j)=A(i,j)*f for i=1,...,nrows and j=1,...,ncols, where A is *this
    // this operator invokes vector_scalar_operation()
    Matrix operator*(Real f) const;
    //! Return <em>B</em> with <em>B</em>(<em>i,j</em>)=<em>A</em>(<em>i,j</em>)/<em>f</em> for <em>i=1,...,nrows</em> and <em>j=1,...,ncols</em>, where <em>A</em> is <em>*this</em>.
    /*! \remark  This operator invokes vector_scalar_operation().
     */
    // return B with B(i,j)=A(i,j)/f for i=1,...,nrows and j=1,...,ncols, where A is *this
    // this operator invokes vector_scalar_operation()
    Matrix operator/(Real f) const;
    //! Return a <b>const</b> reference to <em>A</em> after computing <em>A</em>(<em>i,j</em>)+=<em>f</em> for <em>i=1,...,nrows</em> and <em>j=1,...,ncols</em>, where <em>A</em> is <em>*this</em>.
    /*! \remark  This operator invokes vector_scalar_operation().
     */
    // return a const reference to A after computing A(i,j)+=f for i=1,...,nrows and j=1,...,ncols, where A is *this
    // this operator invokes vector_scalar_operation()
    const Matrix& operator+=(Real f);
    //! Return a <b>const</b> reference to <em>A</em> after computing <em>A</em>(<em>i,j</em>)-=<em>f</em> for <em>i=1,...,nrows</em> and <em>j=1,...,ncols</em>, where <em>A</em> is <em>*this</em>.
    /*! \remark  This operator (eventually) invokes vector_scalar_operation().
     */
    // return a const reference to A after computing A(i,j)-=f for i=1,...,nrows and j=1,...,ncols, where A is *this
    // this operator (eventually) invokes vector_scalar_operation()
    const Matrix& operator-=(Real f);
    //! Return a <b>const</b> reference to <em>A</em> after computing <em>A</em>(<em>i,j</em>)*=<em>f</em> for <em>i=1,...,nrows</em> and <em>j=1,...,ncols</em>, where <em>A</em> is <em>*this</em>.
    /*! \remark  This operator invokes vector_scalar_operation().
     */
    // return a const reference t1o A after computing A(i,j)*=f for i=1,...,nrows and j=1,...,ncols, where A is *this
    // this operator invokes vector_scalar_operation()
    const Matrix& operator*=(Real f);
    //! Return a <b>const</b> reference to <em>A</em> after computing <em>A</em>(<em>i,j</em>)/=<em>f</em> for <em>i=1,...,nrows</em> and <em>j=1,...,ncols</em>, where <em>A</em> is <em>*this</em>.
    /*! \remark  This operator invokes vector_scalar_operation().
     */
    // return a const reference to A after computing A(i,j)/=f for i=1,...,nrows and j=1,...,ncols, where A is *this
    // this operator invokes vector_scalar_operation()
    const Matrix& operator/=(Real f);


    // matrix - matrix addition/subtraction
    //! Return <em>A+B</em>, where <em>A</em> is <em>*this</em>.
    /*! \remark  This operator invokes elementwise_addition().
     */
    // return A+B, where A is *this
    // this operator invokes elementwise_addition()
    Matrix operator+(const Matrix &B) const;
    //! Return <em>A-B</em>, where <em>A</em> is <em>*this</em>.
    /*! \remark  This operator invokes elementwise_addition().
     */
    // return A-B, where A is *this
    // this operator invokes elementwise_addition()
    Matrix operator-(const Matrix &B) const;
    //! Return a <b>const</b> reference to <em>A</em> after computing <em>A=A+B</em>, where <em>A</em> is <em>*this</em>.
    /*! \remark  This operator invokes elementwise_addition().
     */
    // return a const reference to A after computing A = A+B, where A is *this
    // this operator invokes elementwise_addition()
    const Matrix& operator+=(const Matrix &B);
    //! Return a <b>const</b> reference to <em>A</em> after computing <em>A=A-B</em>, where <em>A</em> is <em>*this</em>.
    /*! \remark  This operator invokes elementwise_addition().
     */
    // return a const reference to A after computing A = A-B, where A is *this
    // this operator invokes elementwise_addition()
    const Matrix& operator-=(const Matrix &B);


    // matrix - symmetric matrix addition/subtraction
    //! Return <em>A+B</em>, where <em>A</em> is <em>*this</em>.
    /*! \remark  This operator invokes matrix_symmatrix_addition().
     */
    // return A+B, where A is *this
    // this operator invokes matrix_symmatrix_addition()
    Matrix operator+(const SymmetricMatrix &B) const;
    //! Return <em>A-B</em>, where <em>A</em> is <em>*this</em>.
    /*! \remark  This operator invokes matrix_symmatrix_addition().
     */
    // return A-B, where A is *this
    // this operator invokes matrix_symmatrix_addition()
    Matrix operator-(const SymmetricMatrix &B) const;
    //! Return a <b>const</b> reference to <em>A</em> after computing <em>A=A+B</em>, where <em>A</em> is <em>*this</em>.
    /*! \remark  This operator invokes matrix_symmatrix_addition().
     */
    // return a const reference to A after computing A=A+B, where A is *this
    // this operator invokes matrix_symmatrix_addition()
    const Matrix& operator+=(const SymmetricMatrix &B);
    //! Return a <b>const</b> reference to <em>A</em> after computing <em>A=A-B</em>, where <em>A</em> is <em>*this</em>.
    /*! \remark  This operator invokes matrix_symmatrix_addition().
     */
    // return a const reference to A after computing A=A-B, where A is *this
    // this operator invokes matrix_symmatrix_addition()
    const Matrix& operator-=(const SymmetricMatrix &B);


    // matrix - sparse matrix addition/subtraction
    //! Return <em>A+B</em>, where <em>A</em> is <em>*this</em>.
    /*! \remark  This operator invokes matrix_spmatrix_addition().
     */
    // return A+B, where A is *this
    // this operator invokes matrix_spmatrix_addition()
    Matrix operator+(const SparseMatrix &B) const;
    //! Return <em>A-B</em>, where <em>A</em> is <em>*this</em>.
    /*! \remark  This operator invokes matrix_spmatrix_addition().
     */
    // return A-B, where A is *this
    // this operator invokes matrix_spmatrix_addition()
    Matrix operator-(const SparseMatrix &B) const;
    //! Return a <b>const</b> reference to <em>A</em> after computing <em>A=A+B</em>, where <em>A</em> is <em>*this</em>.
    /*! \remark  This operator invokes matrix_spmatrix_addition().
     */
    // return a const reference to A after computing A=A+B, where A is *this
    // this operator invokes matrix_spmatrix_addition()
    const Matrix& operator+=(const SparseMatrix &B);
    //! Return a <b>const</b> reference to <em>A</em> after computing <em>A=A-B</em>, where <em>A</em> is <em>*this</em>.
    /*! \remark  This operator invokes matrix_spmatrix_addition().
     */
    // return a const reference to A after computing A=A-B, where A is *this
    // this operator invokes matrix_spmatrix_addition()
    const Matrix& operator-=(const SparseMatrix &B);


    // matrix - vector multiplication
    //! Return <em>A*v</em>, where <em>A</em> is <em>*this</em>.
    /*! \remark  This operator invokes matrix_vector_multiplication().
     */
    // return A*v, where A is *this
    // this operator invokes matrix_vector_multiplication()
    Vector operator*(const Vector &v) const;

    // matrix - matrix multiplication
    //! Return <em>A*B</em>, where <em>A</em> is <em>*this</em>.
    /*! \remark  This operator invokes matrix_matrix_multiplication().
     */
    // return A*B, where A is *this
    // this operator invokes matrix_matrix_multiplication()
    Matrix operator*(const Matrix &B) const;
    //! Return a <b>const</b> reference to <em>A</em> after computing <em>A=A*B</em>, where <em>A</em> is <em>*this</em>.
    /*! \remark  This operator (eventually) invokes matrix_matrix_multiplication().
     */
    // return a const reference to A after computing A = A*B, where A is *this
    // this operator (eventually) invokes matrix_matrix_multiplication()
    const Matrix& operator*=(const Matrix &B);


    // matrix - symmetric matrix multiplication
    //! Return <em>A*B</em>, where <em>A</em> is <em>*this</em>.
    /*! \remark  This operator invokes matrix_symmatrix_multiplication().
     */
    // return A*B, where A is *this
    // this operator invokes matrix_symmatrix_multiplication()
    Matrix operator*(const SymmetricMatrix &B) const;
    //! Return a <b>const</b> reference to <em>A</em> after computing <em>A=A*B</em>, where <em>A</em> is <em>*this</em>.
    /*! \remark  This operator invokes matrix_symmatrix_multiplication().
     */
    // return a const reference to A after computing A = A*B, where A is *this
    // this operator (eventually) invokes matrix_symmatrix_multiplication()
    const Matrix& operator*=(const SymmetricMatrix &B);


    // matrix - sparse matrix multiplication
    //! Return <em>A*B</em>, where <em>A</em> is <em>*this</em>.
    /*! \remark  This operator invokes matrix_spmatrix_multiplication().
     */
    // return A*B, where A is *this
    // this operator invokes matrix_spmatrix_multiplication()
    Matrix operator*(const SparseMatrix &B) const;
    //! Return a <b>const</b> reference to <em>A</em> after computing <em>A=A*B</em>, where <em>A</em> is <em>*this</em>.
    /*! \remark  This operator invokes matrix_spmatrix_multiplication().
     */
    // return a const reference to A after computing A = A*B, where A is *this
    // this operator (eventually) invokes matrix_spmatrix_multiplication()
    const Matrix& operator*=(const SparseMatrix &B);


    //! Permute columns of <em>*this</em>.
    /*! Let <em>A</em> be <em>*this</em> and <em>B</em> is the permuted matrix.
     *  <ul>
     *  <li>  If <em>indexFrom1</em>==<b>true</b>,  <em>B</em>(<em>i,perm</em>[<em>j-1</em>])<em> = A</em>(<em>i,j</em>)            for <em>i=1,...,nrows</em> and <em>j=1,...,ncols</em>.
     *  <li>  If <em>indexFrom1</em>==<b>false</b>, <em>B</em>(<em>i,perm</em>[<em>j-1</em>]+<em>1</em>)<em> = A</em>(<em>i,j</em>) for <em>i=1,...,nrows</em> and <em>j=1,...,ncols</em>.
     *  </ul>
     *
     *  \param  perm[]      is the permutation array.
     *  \param  indexFrom1  tells whether the indices start from <em>1</em> or <em>0</em>.
     *
     *  The indices <em>perm</em>[<em>0,...,ncols-1</em>] range over <em>0,1,...,ncols-1</em> (if <em>indexFrom1</em>==<b>false</b>) or over <em>1,2,...,ncols</em> (if <em>indexFrom1</em>==<b>true</b>).
     *
     *  \return The permuted matrix <em>B</em>.
     *
     *  \remark  This member function invokes permute_matrix_columns().
     *  \remark  This member function resembles <em>B</em>(:,<em>perm</em>)=<em>A</em> in <b>OCTAVE</b>/<b>MATLAB</b>.
     */
    // permute columns of A, such that
    //    B(i,perm[j-1])   = A(i,j) for i=1,...,nrows and j=1,...,ncols (if indexFrom1==true), or
    //    B(i,perm[j-1]+1) = A(i,j) for i=1,...,nrows and j=1,...,ncols (if indexFrom1==false),
    // where A is *this and B is the return matrix
    // the indices perm[0,...,ncols-1] range over 0,1,...,ncols-1 (if indexFrom1==false) or over 1,2,...,ncols (if indexFrom1==true)
    // this member function invokes permute_matrix_columns()
    // this member function resembles B(:,perm) = A in OCTAVE/MATLAB
    Matrix permuteColumns(const mkIndex *perm, bool indexFrom1=true) const;

    //! Permute rows of <em>*this</em>.
    /*! Let <em>A</em> be <em>*this</em> and <em>B</em> is the permuted matrix.
     *  <ul>
     *  <li>  If <em>indexFrom1</em>==<b>true</b>,  <em>B</em>(<em>perm</em>[<em>i-1</em>]<em>,j</em>)<em>   = A</em>(<em>i,j</em>) for <em>i=1,...,nrows</em> and <em>j=1,...,ncols</em>.
     *  <li>  If <em>indexFrom1</em>==<b>false</b>, <em>B</em>(<em>perm</em>[<em>i-1</em>]+<em>1,j</em>)<em> = A</em>(<em>i,j</em>) for <em>i=1,...,nrows</em> and <em>j=1,...,ncols</em>.
     *  </ul>
     *  \param  perm[]      is the permutation array.
     *  \param  indexFrom1  tells whether the indices start from <em>1</em> or <em>0</em>.
     *
     *  The indices <em>perm</em>[<em>0,...,nrows-1</em>] range over <em>0,1,...,nrows-1</em> (if <em>indexFrom1</em>==<b>false</b>) or over <em>1,2,...,nrows</em> (if <em>indexFrom1</em>==<b>true</b>).
     *
     *  \return  The permuted matrix <em>B</em>.
     *
     *  \remark  This member function invokes permute_matrix_rows().
     *  \remark  This member function resembles <em>B</em>(<em>perm</em>,:)=<em>A</em> in <b>OCTAVE</b>/<b>MATLAB</b>.
     */
    // permute rows of A, such that
    //    B(perm[i-1],j)   = A(i,j) for i=1,...,nrows and j=1,...,ncols (if indexFrom1==true), or
    //    B(perm[i-1]+1,j) = A(i,j) for i=1,...,nrows and j=1,...,ncols (if indexFrom1==false),
    // where A is *this and B is the return matrix
    // the indices perm[0,...,nrows-1] range over 0,1,...,nrows-1 (if indexFrom1==false) or over 1,2,...,nrows (if indexFrom1==true)
    // this member function invokes permute_matrix_rows()
    // this member function resembles B(perm,:) = A in OCTAVE/MATLAB
    Matrix permuteRows(const mkIndex *perm, bool indexFrom1=true) const;

    //! Permute rows and columns of <em>*this</em>.
    /*! Let <em>A</em> be <em>*this</em> and <em>B</em> is the permuted matrix.
     *  <ul>
     *  <li>  If <em>indexFrom1</em>==<b>true</b>, <em>B</em>(<em>rperm</em>[<em>i-1</em>],<em>cperm</em>[<em>j-1</em>])               = <em>A</em>(<em>i,j</em>) for <em>i=1,...,nrows</em> and <em>j=1,...,ncols</em>.
     *  <li>  If <em>indexFrom1</em>==<b>false</b>, <em>B</em>(<em>rperm</em>[<em>i-1</em>]+<em>1,cperm</em>[<em>j-1</em>]+<em>1</em>) = <em>A</em>(<em>i,j</em>) for <em>i=1,...,nrows</em> and <em>j=1,...,ncols</em>.
     *  </ul>
     *  \param  rperm[]     is the row permutation array.
     *  \param  cperm[]     is the column permutation array.
     *  \param  indexFrom1  tells whether the indices start from <em>1</em> or <em>0</em>.
     *
     *  The indices <em>rperm</em>[<em>0,...,nrows-1</em>] range over <em>0,1,...,nrows-1</em> (if <em>indexFrom1</em>==<b>false</b>) or over <em>1,2,...,nrows</em> (if <em>indexFrom1</em>==<b>true</b>). \n
     *  The indices <em>cperm</em>[<em>0,...,ncols-1</em>] range over <em>0,1,...,ncols-1</em> (if <em>indexFrom1</em>==<b>false</b>) or over <em>1,2,...,ncols</em> (if <em>indexFrom1</em>==<b>true</b>).
     *
     *  \remark  If <em>cperm</em>==<tt>NULL</tt>, then symmetric permutation will be applied (i.e., <em>rperm</em>[] will also be used as <em>cperm</em>[]).
     *  \remark  This member function (eventually) invokes permute_matrix_rows() and permute_matrix_columns().
     *  \remark  This member function resembles <em>B</em>(<em>rperm,cperm</em>)=<em>A</em> in <b>OCTAVE</b>/<b>MATLAB</b>.
     */
    // permute rows and columns of A, such that
    //    B(rperm[i-1],cperm[j-1])     = A(i,j) for i=1,...,nrows and j=1,...,ncols (if indexFrom1==true), or
    //    B(rperm[i-1]+1,cperm[j-1]+1) = A(i,j) for i=1,...,nrows and j=1,...,ncols (if indexFrom1==false),
    // where A is *this and B is the return matrix
    // if cperm==NULL, then symmetric permutation will be applied (i.e., rperm[] will also be used as cperm[])
    // the indices rperm[0,...,nrows-1] range over 0,1,...,nrows-1 (if indexFrom1==false) or over 1,2,...,nrows (if indexFrom1==true)
    // the indices cperm[0,...,ncols-1] range over 0,1,...,ncols-1 (if indexFrom1==false) or over 1,2,...,ncols (if indexFrom1==true)
    // this member function (eventually) invokes permute_matrix_rows() and permute_matrix_columns()
    // this member function resembles B(rperm,cperm) = A in OCTAVE/MATLAB
    Matrix permuteRowsAndColumns(const mkIndex *rperm, const mkIndex *cperm=NULL, bool indexFrom1=true) const;


    #ifdef USE_MEX
    // convert *mxArray to Matrix
    static Matrix mxArray2Matrix(const mxArray *mA);
    // convert Matrix to *mxArray
    static mxArray *Matrix2mxArray(const Matrix &A);
    #endif
};


// matrix scalar operations
//! Return <em>f*A</em>.
/*! This operator (eventually) invokes vector_scalar_operation().
 */
// return f*A
// this operator (eventually) invokes vector_scalar_operation()
Matrix operator*(Real h, const Matrix &A);

//! Return <em>B</em> with <em>B</em>(<em>i,j</em>)=<em>f+A</em>(<em>i,j</em>) for <em>i=1,...,A.nrows</em> and <em>j=1,...,A.ncols</em>.
/*! This operator (eventually) invokes vector_scalar_operation().
 */
// return B with B(i,j)=f+A(i,j) for i=1,...,A.nrows and j=1,...,A.ncols
// this operator (eventually) invokes vector_scalar_operation()
Matrix operator+(Real h, const Matrix &A);

//! Return <em>B</em> with <em>B</em>(<em>i,j</em>)=<em>f-A</em>(<em>i,j</em>) for <em>i=1,...,A.nrows</em> and <em>j=1,...,A.ncols</em>.
/*! This operator (eventually) invokes vector_scalar_operation().
 */
// return B with B(i,j)=f-A(i,j) for i=1,...,A.nrows and j=1,...,A.ncols
// this operator (eventually) invokes vector_scalar_operation()
Matrix operator-(Real h, const Matrix &A);


//! Display matrix <em>A</em> with operator <<.
/*! In other words, print out matrix <em>A</em> to <b>ostream</b> <em>s</em>.
 *  \remark  This routine invokes print_matrix().
 */
// display matrix A with operator <<
// in other words, print out matrix A to ostream s
// this routine invokes print_matrix()
std::ostream& operator<<(std::ostream &s, const Matrix &A);


#ifdef USE_NAMESPACE
}  // end of namespace MATKIT
#endif

#endif  // end of #ifdef MATRIX_H
