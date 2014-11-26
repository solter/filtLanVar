#ifndef SPARSE_MATRIX_H
#define SPARSE_MATRIX_H
/*! \file spmatrix.h
 *  The class SparseMatrix defines a sparse matrix
 *  stored in compressed sparse column (CSC) format in
 *  SparseMatrix::store[], SparseMatrix::rowIndex[], and SparseMatrix::colPointer[].
 *
 *  <ul>
 *  <li>  The number of nonzero elements in <em>j</em>th column is SparseMatrix::colPointer[<em>j</em>] - SparseMatrix::colPointer[<em>j-1</em>].
 *  <li>  The nonzero elements in <em>j</em>th column are stored in SparseMatrix::store[SparseMatrix::colPointer[<em>j-1</em>],...,SparseMatrix::colPointer[<em>j</em>]-1].
 *  <li>  The corresponding row indices are stored in SparseMatrix::rowIndex[SparseMatrix::colPointer[<em>j-1</em>],...,SparseMatrix::colPointer[<em>j</em>]-1]. \n
 *        (ps. The index of <em>i</em>th row is <em>i-1</em>.)
 *  </ul>
 *  In the above we have assumed SparseMatrix::colPointer[<em>0</em>]==<em>0</em> for simplicity.
 */

#include <iostream>  // for ostream, etc.
#include "matkitdef.h"


#ifdef USE_NAMESPACE
namespace MATKIT {
#endif


class Vector;
class Matrix;
class SymmetricMatrix;

/*! \brief Class SparseMatrix defines a generic matrix and the arithmetic operations on it,
 *         and also providing some static functions for sparse matrix generation.
 */
class SparseMatrix {
protected:
    // number of rows and number of columns of *this
    //! Number of rows of <em>*this</em>.
    // number of rows of *this
    mkIndex nrows;
    //! Number of columns of <em>*this</em>.
    // number of columns of *this
    mkIndex ncols;

    // CSC format
    //! A storage for storing values of nonzero elements of <em>*this</em>.
    // a storage for storing values of nonzero elements of *this
    Real *store;
    //! An array for row indices of <em>*this</em>.
    // an array for row indices of *this
    mkIndex *rowIndex;
    //! An array for column pointers of <em>*this</em>.
    // an array for column pointers of *this
    mkIndex *colPointer;
    //! Maximum number of nonzero elements allowed in the allocated memory.
    // maximum number of nonzero elements allowed in the allocated memory
    mkIndex maxnnz;

public:
    //! When the allocated memory is insufficient, this is the factor by which the memory is expanded (default <em>1.2</em>).
    /*! \remark  This parameter is used by operator()(<b>mkIndex</b>,<b>mkIndex</b>).
     *           To be precise, when setting <em>A</em>(<em>i,j</em>)=<em>f</em> with <em>A</em> being <em>*this</em>,
     *           an additional space is needed for storing <em>A</em>(<em>i,j</em>) if <em>A</em>(<em>i,j</em>) is previously zero.
     *           When memory is insufficient, it requires an expansion.
     *  \see  SparseMatrix::setMemoryExpansionFactor().
     */
    // when the allocated memory is insufficient, this is the factor by which the memory is expanded (default 1.2)
    // this parameter is used by operator()(mkIndex,mkIndex)
    // to be precise, when setting A(i,j)=f with A being *this,
    // an additional space is needed for storing A(i,j) if A(i,j) is previously zero
    // when memory is insufficient, it requires an expansion
    // see also SparseMatrix::setMemoryExpansionFactor()
    static Real memoryExpansionFactor;

    //! Set the memory expansion factor.
    /*! The memory expansion factor determines how much memory should be reallocated
     *  if the allocated memory is insufficient, which can happen when SparseMatrix::operator()(<b>mkIndex</b>,<b>mkIndex</b>) is called.
     *  \see SparseMatrix::memoryExpansionFactor.
     */
    // set memory expansion factor
    // the memory expansion factor determines how much memory should be reallocated
    // if the allocated memory is insufficient, which can happen when SparseMatrix::operator()(mkIndex,mkIndex) is called
    // the default value of memoryExpansionFactor is 1.2, shown in spmatrix.cpp
    // see also SparseMatrix::memoryExpansionFactor
    static void setMemoryExpansionFactor(Real newExpansionFactor);

    //! Read a sparse matrix from a matrix-market file <em>fileName</em>[] and return the result.
    /*! \param  fileName[]  is the file name of the matrix-market file.
     *  \param  sort        tells whether to sort elements in each column with respect to the row indices. \n
     *                      Usually a matrix-market file has elements sorted, in which case <em>sort</em>==<b>true</b> is unnecessary.
     *  \remark  This member function invokes matrix_market_spmatrix_read().
     */
    // read a sparse matrix from a matrix-market file and return the result
    // fileName[] is the file name of the matrix-market file
    // sort       tells whether to sort elements in each column with respect to the row indices
    //            usually a matrix-market file has elements sorted, in which case sort==true is unnecessary
    // this member function invokes matrix_market_spmatrix_read()
    static SparseMatrix mmread(const char fileName[], bool sort=true);

    //! Write a sparse matrix <em>A</em> to a matrix-market file <em>fileName</em>[], where <em>A</em> is <em>*this</em>.
    /*! \param  fileName[]   is a character string storing the output matrix-market file name.
     *  \param  patternOnly  tells whether only the sparsity pattern is considered. \n
     *                       If <em>patternOnly</em>==<b>true</b>, then only the sparsity pattern is to be written to the file <em>fileName</em>[], and <em>A.store</em>[] will not be used.
     *  \param  comment[]    is a character string as a comment. \n
     *                       If \f$\neq\f$<tt>NULL</tt>, then the comment will be printed in the output matrix-market file.
     *
     *  \return  An integer which indicates the status:
     *  <ul>
     *  <li>  <em>0</em>:  a successful write.
     *  <li>  <em>1</em>:  file open error.
     *  <li>  <em>2</em>:  file close error.
     *  <li>  A negative integer:  a number returned by <b>fprintf</b>() which signals an error.
     *  </ul>
     *
     *  \remark  This member function invokes matrix_market_spmatrix_write().
     */
    // write a sparse matrix A to a matrix-market file fileName[], where A is *this
    // if patternOnly==true, then only the sparsity pattern is to be written to the file fileName[], and
    //     A.store[] will not be used
    // comment[] is a character string which will be printed in the output file fileName[] (if comment!=NULL)
    // the return value:
    //    0: a successful read
    //    1: file open error
    //    2: file close error
    //    a negative integer: a number returned by fprintf() which signals an error
    // this member function invokes matrix_market_spmatrix_write()
    int mmwrite(const char fileName[], bool patternOnly=false, const char comment[]=NULL) const {
        return matrix_market_spmatrix_write(nrows, ncols, store, rowIndex, colPointer, fileName, patternOnly, comment);
    }


    // some static functions
    //! Return an <em>n</em>-by-<em>n</em> (sparse) matrix of zeros.
    // return an n-by-n (sparse) matrix of zeros
    static SparseMatrix zeros(mkIndex n) { return SparseMatrix(n, n); }
    //! Return an <em>nr</em>-by-<em>nc</em> (sparse) matrix of zeros.
    // return an nr-by-nc (sparse) matrix of zeros
    static SparseMatrix zeros(mkIndex nr, mkIndex nc) { return SparseMatrix(nr, nc); }
    //! Return an <em>n</em>-by-<em>n</em> (sparse) matrix with ones on the main diagonal and zeros elsewhere.
    /*! \remark  This static function resembles <b>speye</b>(<em>n</em>) in <b>OCTAVE</b>/<b>MATLAB</b>.
     */
    // return an n-by-n (sparse) matrix with ones on the main diagonal and zeros elsewhere
    // this static function resembles speye(n) in OCTAVE/MATLAB
    static SparseMatrix eye(mkIndex n) { return eye(n,n); }
    //! Return an <em>nr</em>-by-<em>nc</em> (sparse) matrix with ones on the main diagonal and zeros elsewhere.
    /*! \remark  This static function resembles <b>speye</b>(<em>nr</em>,<em>nc</em>) in <b>OCTAVE</b>/<b>MATLAB</b>.
     */
    // return an nr-by-nc (sparse) matrix with ones on the main diagonal and zeros elsewhere
    // this static function resembles speye(nr,nc) in OCTAVE/MATLAB
    static SparseMatrix eye(mkIndex nr, mkIndex nc);

    //! Return a tridiagonal matrix <em>A</em> as a sparse matrix with three diagonals <em>ldiag</em>[], <em>diag</em>[], and <em>udiag</em>[].
    /*! The lower subdiagonal, the main diagonal, and the upper subdiagonal of <em>A</em> are formed by
     *  <em>ldiag</em>[<em>0,...,nrc-2</em>], <em>diag</em>[<em>0,...,nrc-1</em>], and <em>udiag</em>[<em>0,...,nrc-2</em>], respectively.
     *  \param  nrc     is the dimension of <em>A</em>.
     *  \param  ldiag[] stores the subdiagonal elements in the lower triangular part of <em>A</em>. \n
     *          <em>A</em>(<em>i+1,i</em>) = <em>ldiag</em>[<em>i-1</em>] for <em>i=1,...,nrc-1</em>.
     *  \param  diag[]  stores the elements of <em>A</em> on the main diagonal. \n
     *          <em>A</em>(<em>i,i</em>)   = <em> diag</em>[<em>i-1</em>] for <em>i=1,...,nrc</em>.
     *  \param  udiag[] stores the subdiagonal elements in the upper triangular part of <em>A</em>. \n
     *          <em>A</em>(<em>i,i+1</em>) = <em>udiag</em>[<em>i</em>]   for <em>i=1,...,nrc-1</em>.
     *  \return  The sparse matrix <em>A</em>.
     *  \remark  This static function invokes tridiagonal_to_spmatrix().
     */
    // return a tridiagonal matrix A as a sparse matrix with three diagonals ldiag[], diag[], and udiag[]
    // the lower subdiagonal, the main diagonal, and the upper subdiagonal of A are formed by
    // ldiag[0,...,nrc-2], diag[0,...,nrc-1], and udiag[0,...,nrc-2], respectively
    // to be specific, A(i+1,i) = ldiag[i-1] for i=1,...,nrc-1,
    //                 A(i,i  ) =  diag[i-1] for i=1,...,nrc,
    //             and A(i,i+1) = udiag[i]   for i=1,...,nrc-1
    // this static function invokes tridiagonal_to_spmatrix()
    static SparseMatrix tridiagonalMatrix(mkIndex nrc, const Real *ldiag, const Real *diag, const Real *udiag);


    // constructors
    //! A constructor for an empty (sparse) matrix.
    // a constructor for an empty (sparse) matrix
    SparseMatrix();
    //! A constructor for an <em>nr</em>-by-<em>nc</em> (sparse) matrix of zeros.
    // a constructor for an nr-by-nc (sparse) matrix of zeros
    SparseMatrix(mkIndex nr, mkIndex nc);
    //! A constructor for an nr-by-nc (sparse) matrix stored in CSC format in <em>store0</em>[], <em>rowIndex0</em>[], <em>colPointer0</em>[].
    /*! \remark  It performs a shallow copy of <em>store0</em>, <em>rowIndex0</em>, and <em>colPointer0</em>, and the memory will be freed by the destructor.
     *  \remark  If <em>maxnnz0</em> is greater than the number of nonzero elements, then it is assumed that
     *           memory of size <em>maxnnz0</em> is allocated with the address pointed to by <em>store0</em>, and so is it for <em>rowIndex0</em>.
     */
    // a constructor for an nr-by-nc (sparse) matrix stored in CSC format in store0[], rowIndex0[], colPointer0[]
    // note that it performs a shallow copy of store0, rowIndex0, and colPointer0, and the memory will be freed by the destructor
    // if maxnnz0 is greater than the number of nonzero elements, then it is assumed that
    // memory of size maxnnz0 is allocated with the address pointed to by store0, and so is it for rowIndex0
    SparseMatrix(mkIndex nr, mkIndex nc, Real *store0, mkIndex *rowIndex0, mkIndex *colPointer0, mkIndex maxnzz=0);
    //! A copy constructor.
    /*! \remark  This constructor invokes memory_xcopy() for copying nonzero elements.
     */
    // a copy constructor
    // this constructor invokes memory_xcopy() for copying nonzero elements
    SparseMatrix(const SparseMatrix &B);

    //! A destructor. Memory allocated, with the addresses pointed to by data members SparseMatrix::store, SparseMatrix::rowIndex, SparseMatrix::colPointer, will be freed.
    // a destructor; memory allocated, with the addresses pointed to by data members SparseMatrix::store, SparseMatrix::rowIndex, SparseMatrix::colPointer, will be freed
    virtual ~SparseMatrix();


    // get the information of *this
    //! The number of rows of <em>*this</em>.
    // the number of rows of *this
    mkIndex Nrows() const { return nrows; }
    //! The number of columns of <em>*this</em>.
    // the number of columns of *this
    mkIndex Ncols() const { return ncols; }
    //! A <b>const</b> pointer to the storage for values of non-zero elements of <em>*this</em>.
    // a const pointer to the storage for values of non-zero elements of *this
    const Real *Store() const { return store; }
    //! A <b>const</b> pointer to the array of row indices of <em>*this</em>.
    // a const pointer to the array of row indices of *this
    const mkIndex *RowIndex() const { return rowIndex; }
    //! A <b>const</b> pointer to the array of column pointers of <em>*this</em>.
    // a const pointer to the array of column pointers of *this
    const mkIndex *ColPointer() const { return colPointer; }
    //! Number of nonzero elements.
    // number of nonzero elements
    mkIndex Nnz() const { return colPointer[ncols]-colPointer[0]; }
    //! Maximum number of nonzero elements allowed in the current allocated memory.
    // maximum number of nonzero elements allowed in the current allocated memory
    mkIndex Maxnnz() const { return maxnnz; }


    //! Convert <em>*this</em> to a general matrix.
    /*! \remark  This static function invokes sparse_to_full().
     *  \remark  This static function resembles <b>full</b>(<em>A</em>) in <b>OCTAVE</b>/<b>MATLAB</b>, where <em>A</em> is <em>*this</em>.
     */
    // convert *this to a general matrix
    // this static function invokes sparse_to_full()
    // this static function resembles full(A) in OCTAVE/MATLAB, where A is *this
    Matrix toFull() const;

    //! Set <em>*this</em> as a sparse matrix of size <em>nr</em>-by-<em>nc</em> stored in CSC format in <em>store0</em>[], <em>rowIndex0</em>[], <em>colPointer0</em>[].
    /*! \remark  It performs a shallow copy of <em>store0</em>, <em>rowIndex0</em>, and <em>colPointer0</em>, and the memory will be freed by the destructor.
     *  \remark  If <em>maxnnz0</em> is greater than the number of nonzero elements, then it is assumed that
     *           memory of size <em>maxnnz0</em> is allocated with the address pointed to by <em>store0</em>, and so is it for <em>rowIndex0</em>.
     */
    // set *this as a sparse matrix of size nr-by-nc stored in CSC format in store0[], rowIndex0[], colPointer0[]
    // note that it performs a shallow copy of store0, rowIndex0, and colPointer0, and the memory will be freed by the destructor
    // if maxnnz0 is greater than the number of nonzero elements, then it is assumed that
    // memory of size maxnnz0 is allocated with the address pointed to by store0, and so is it for rowIndex0
    void set(mkIndex nr, mkIndex nc, Real *store0, mkIndex *rowIndex0, mkIndex *colPointer0, mkIndex maxnnz0=0);


    //! A submatrix formed by columns <em>j1,...,j2</em> of <em>*this</em>.
    /*! \remark  This member function invokes submatrix_of_matrix().
     *  \remark  This member function resembles <em>A</em>(:,<em>j1:j2</em>) in <b>OCTAVE</b>/<b>MATLAB</b>, where <em>A</em> is <em>*this</em>.
     *  \see  SparseMatrix::rows(), SparseMatrix::subMatrix().
     */
    // a submatrix formed by columns j1,...,j2 of *this
    // this member function invokes submatrix_of_spmatrix()
    // this member function resembles A(:,j1:j2) in OCTAVE/MATLAB, where A is *this
    SparseMatrix columns(mkIndex j1, mkIndex j2) const;
    //! A submatrix formed by rows <em>i1,...,i2</em> of <em>*this</em>.
    /*! \remark  This member function invokes submatrix_of_spmatrix().
     *  \remark  This member function resembles <em>A</em>(<em>i1:i2</em>,:) in <b>OCTAVE</b>/<b>MATLAB</b>, where <em>A</em> is <em>*this</em>.
     *  \see  SparseMatrix::columns(), SparseMatrix::subMatrix().
     */
    // a submatrix formed by rows i1,...,i2 of *this
    // this member function invokes submatrix_of_spmatrix()
    // this member function resembles A(i1:i2,:) in OCTAVE/MATLAB, where A is *this
    SparseMatrix rows(mkIndex i1, mkIndex i2) const;
    //! A submatrix formed by elements in rows <em>i1,...,i2</em> and columns <em>j1,...,j2</em> of <em>*this</em>.
    /*! \remark  This member function invokes submatrix_of_spmatrix().
     *  \remark  This member function resembles <em>A</em>(<em>i1:i2,j1:j2</em>) in <b>OCTAVE</b>/<b>MATLAB</b>, where <em>A</em> is <em>*this</em>.
     *  \see  SparseMatrix::rows(), SparseMatrix::columns().
     */
    // a submatrix formed by elements in rows i1,...,i2 and columns j1,...,j2 of *this
    // this member function invokes submatrix_of_matrix()
    // this member function resembles A(i1:i2,j1:j2) in OCTAVE/MATLAB, where A is *this
    SparseMatrix subMatrix(mkIndex i1, mkIndex i2, mkIndex j1, mkIndex j2) const;


    //! A vector formed by the elements in the <em>k</em>th diagonal of <em>*this</em> in order.
    /*! \remark  This member function invokes matrix_diag().
     *  \remark  This member function resembles <b>diag</b>(<em>A</em>,<em>k</em>) in <b>OCTAVE</b>/<b>MATLAB</b>.
     */
    // a vector formed by the elements in the k-th diagonal of *this in order
    // this member function invokes spmatrix_diag()
    // this member function resembles diag(A,k) in OCTAVE/MATLAB
    Vector diag(mkSignedIndex k=0) const;
    //! The transpose of <em>*this</em>.
    /*! \remark  This member function invokes spmatrix_transpose().
     *  \remark  This member function resembles <em>A</em>' or equivalently <b>transpose</b>(<em>A</em>) in <b>OCTAVE</b>/<b>MATLAB</b>,
     *           where <em>A</em> is <em>*this</em>.
     */
    // the transpose of *this
    // this member function invokes spmatrix_transpose()
    // this member function resembles A' or equivalently transpose(A) in OCTAVE/MATLAB, where A is *this
    SparseMatrix transpose() const;


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
    /*! \remark  If <em>A</em>(<em>i,j</em>)==<em>0</em>, then a space will be created for it.
     */
    // return a reference to A(i,j), i.e. can be used to write A(i,j), where A is *this
    // if A(i,j)==0, then a space will be created for it
    virtual Real& operator()(mkIndex i, mkIndex j);
    //! Read <em>A</em>(<em>i,j</em>) via a <b>const</b> reference or pointer to <em>A</em>, where <em>A</em> is <em>*this</em>.
    // read A(i,j) via a const reference or pointer to A, where A is *this
    Real operator()(mkIndex i, mkIndex j) const;


    // operators ==, !=
    //! Return <b>true</b> if <em>A==B</em> mathematically, and otherwise return <b>false</b>, where <em>A</em> is <em>*this</em>.
    /*! \remark  This operator invokes is_spmatrix_equal_spmatrix().
     */
    // return true if A == B mathematically, and otherwise return false, where A is *this
    // this operator invokes is_spmatrix_equal_spmatrix()
    bool operator==(const SparseMatrix &B) const {
        if (nrows != B.nrows || ncols != B.ncols)
            return false;
        return is_spmatrix_equal_spmatrix(ncols, Store(), RowIndex(), ColPointer(), B.Store(), B.RowIndex(), B.ColPointer());
    }
    //! Return <b>false</b> if <em>A==B</em> mathematically, and otherwise return <b>true</b>, where <em>A</em> is <em>*this</em>.
    /*! \remark  This operator (eventually) invokes is_spmatrix_equal_spmatrix().
     */
    // return false if A == B mathematically, and otherwise return true, where A is *this
    // this operator (eventually) invokes is_spmatrix_equal_spmatrix()
    bool operator!=(const SparseMatrix &B) const {
        return !(*this==B);
    }


    // operator =
    //! Return a <b>const</b> reference to <em>A</em> after setting <em>A=B</em>, where <em>A</em> is <em>*this</em>.
    /*! This operator invokes memory_xcopy() for copying nonzero elements.
     */
    // return a const reference to A after setting A=B, where A is *this
    // this operator invokes memory_xcopy() for copying nonzero elements
    const SparseMatrix& operator=(const SparseMatrix &B);


    // operator -
    //! Return -<em>A</em>, where <em>A</em> is <em>*this</em>.
    /*! This operator invokes vector_scalar_operation() for copying elements.
     */
    // return -A, where A is *this
    // this operator invokes vector_scalar_operation() for copying elements
    SparseMatrix operator-() const;


    // operator *,/,*=,/= a real number
    //! Return <em>B</em> with <em>B</em>(<em>i,j</em>)=<em>A</em>(<em>i,j</em>)*<em>f</em> for <em>i,j=1,...,nrows</em>, where <em>A</em> is <em>*this</em>.
    /*! \remark  This operator invokes vector_scalar_operation() to copy scaled nonzero elements.
     */
    // return B with B(i,j)=A(i,j)*f for i,j=1,...,nrows, where A is *this
    // this operator invokes vector_scalar_operation() to copy scaled nonzero elements
    SparseMatrix operator*(Real f) const;
    //! Return <em>B</em> with <em>B</em>(<em>i,j</em>)=<em>A</em>(<em>i,j</em>)/<em>f</em> for <em>i,j=1,...,nrows</em>, where <em>A</em> is <em>*this</em>.
    /*! \remark  This operator invokes vector_scalar_operation() to copy scaled nonzero elements.
     */
    // return B with B(i,j)=A(i,j)/f for i,j=1,...,nrows, where A is *this
    // this operator invokes vector_scalar_operation() to copy scaled nonzero elements
    SparseMatrix operator/(Real f) const;
    //! Return a <b>const</b> reference to <em>A</em> after computing <em>A</em>(<em>i,j</em>)*=<em>f</em> for <em>i,j=1,...,nrows</em>, where <em>A</em> is <em>*this</em>.
    /*! \remark  This operator invokes vector_scalar_operation() to scale nonzero elements.
     */
    // return a const reference to A after computing A(i,j)*=f for i=1,...,nrows and j=1,...,ncols, where A is *this
    // this operator invokes vector_scalar_operation() to scale nonzero elements
    const SparseMatrix& operator*=(Real f);
    //! Return a <b>const</b> reference to <em>A</em> after computing <em>A</em>(<em>i,j</em>)/=<em>f</em> for <em>i,j=1,...,nrows</em>, where <em>A</em> is <em>*this</em>.
    /*! \remark  This operator invokes vector_scalar_operation() to scale nonzero elements.
     */
    // return a const reference to A after computing A(i,j)/=f for i=1,...,nrows and j=1,...,ncols, where A is *this
    // this operator invokes vector_scalar_operation() to scale nonzero elements
    const SparseMatrix& operator/=(Real f);


    // operator +,-,+=,-= a sparse matrix
    //! Return <em>A+B</em>, where <em>A</em> is <em>*this</em>.
    /*! \remark  This operator invokes spmatrix_spmatrix_addition().
     */
    // return A+B, where A is *this
    // this member function invokes spmatrix_spmatrix_addition()
    SparseMatrix operator+(const SparseMatrix &B) const;
    //! Return <em>A-B</em>, where <em>A</em> is <em>*this</em>.
    /*! \remark  This operator invokes spmatrix_spmatrix_addition().
     */
    // return A-B, where A is *this
    // this member function invokes spmatrix_spmatrix_addition()
    SparseMatrix operator-(const SparseMatrix &B) const;
    //! Return a <b>const</b> reference to <em>A</em> after computing <em>A=A+B</em>, where <em>A</em> is <em>*this</em>.
    /*! \remark  This operator (eventually) invokes spmatrix_spmatrix_addition().
     */
    // return a const reference to A after computing A = A+B, where A is *this
    // this member function (eventually) invokes spmatrix_spmatrix_addition()
    const SparseMatrix& operator+=(const SparseMatrix &B);
    //! Return a <b>const</b> reference to <em>A</em> after computing <em>A=A-B</em>, where <em>A</em> is <em>*this</em>.
    /*! \remark  This operator (eventually) invokes spmatrix_spmatrix_addition().
     */
    // return a const reference to A after computing A = A-B, where A is *this
    // this member function (eventually) invokes spmatrix_spmatrix_addition()
    const SparseMatrix& operator-=(const SparseMatrix &B);


    // sparse matrix - matrix addition/subtraction
    //! Return <em>A+B</em>, where <em>A</em> is <em>*this</em>.
    /*! \remark  This operator (eventually) invokes matrix_spmatrix_addition().
     */
    // return A+B, where A is *this
    // this member function (eventually) invokes matrix_spmatrix_addition()
    Matrix operator+(const Matrix &B) const;
    //! Return <em>A-B</em>, where <em>A</em> is <em>*this</em>.
    /*! \remark  This operator (eventually) invokes matrix_spmatrix_addition().
     */
    // return A-B, where A is *this
    // this member function (eventually) invokes matrix_spmatrix_addition()
    Matrix operator-(const Matrix &B) const;


    // sparse matrix - vector multiplication
    //! Return <em>A*v</em>, where <em>A</em> is <em>*this</em>.
    /*! \remark  This operator invokes spmatrix_vector_multiplication().
     */
    // return A*v, where A is *this
    // this operator invokes spmatrix_vector_multiplication()
    Vector operator*(const Vector &v) const;


    // sparse matrix - matrix multiplication
    //! Return <em>A*B</em>, where <em>A</em> is <em>*this</em>.
    /*! \remark  This operator invokes spmatrix_matrix_multiplication().
     */
    // return A*B, where A is *this
    // this operator invokes spmatrix_matrix_multiplication()
    Matrix operator*(const Matrix &B) const;

    //! Permute columns of <em>*this</em>.
    /*! Let <em>A</em> be <em>*this</em> and <em>B</em> is the permuted matrix.
     *  <ul>
     *  <li>  If <em>indexFrom1</em>==<b>true</b>,  <em>B</em>(<em>i,perm</em>[<em>j-1</em>])<em> = A</em>(<em>i,j</em>)            for <em>i=1,...,nrows</em> and <em>j=1,...,ncols</em>.
     *  <li>  If <em>indexFrom1</em>==<b>false</b>, <em>B</em>(<em>i,perm</em>[<em>j-1</em>]+<em>1</em>)<em> = A</em>(<em>i,j</em>) for <em>i=1,...,nrows</em> and <em>j=1,...,ncols</em>.
     *  </ul>
     *  \param  perm[] is the permutation array.
     *  \param  indexFrom1 tells whether the indices start from <em>1</em> or <em>0</em>.
     *
     *  The indices <em>perm</em>[<em>0,...,ncols-1</em>] range over <em>0,1,...,ncols-1</em> (if <em>indexFrom1</em>==<b>false</b>) or over <em>1,2,...,ncols</em> (if <em>indexFrom1</em>==<b>true</b>).
     *
     *  \return  The permuted matrix <em>B</em>.
     *  \remark  This member function invokes permute_csc_columns().
     *  \remark  This member function resembles <em>B</em>(:,<em>perm</em>)=<em>A</em> in <b>OCTAVE</b>/<b>MATLAB</b>.
     */
    // permute columns of A, such that
    //    B(i,perm[j-1])   = A(i,j) for i=1,...,nrows and j=1,...,ncols (if indexFrom1==true), or
    //    B(i,perm[j-1]+1) = A(i,j) for i=1,...,nrows and j=1,...,ncols (if indexFrom1==false), where
    //    A is *this and B is the return sparse matrix
    // the indices cperm[0,...,nrows-1] range over 0,1,...,ncols-1 (if indexFrom1==false) or over 1,2,...,ncols (if indexFrom1==true)
    // this member function invokes permute_csc_columns()
    // this member function resembles B(:,perm) = A in OCTAVE/MATLAB
    SparseMatrix permuteColumns(const mkIndex *perm, bool indexFrom1=true) const;

    //! Permute rows of <em>*this</em>.
    /*! Let <em>A</em> be <em>*this</em> and <em>B</em> is the permuted matrix.
     *  <ul>
     *  <li>  If <em>indexFrom1</em>==<b>true</b>,  <em>B</em>(<em>perm</em>[<em>i-1</em>]<em>,j</em>)<em>   = A</em>(<em>i,j</em>) for <em>i=1,...,nrows</em> and <em>j=1,...,ncols</em>.
     *  <li>  If <em>indexFrom1</em>==<b>false</b>, <em>B</em>(<em>perm</em>[<em>i-1</em>]+<em>1,j</em>)<em> = A</em>(<em>i,j</em>) for <em>i=1,...,nrows</em> and <em>j=1,...,ncols</em>.
     *  </ul>
     *  \param  perm[] is the permutation array.
     *  \param  indexFrom1 tells whether the indices start from <em>1</em> or <em>0</em>.
     *
     *  The indices <em>perm</em>[<em>0,...,nrows-1</em>] range over <em>0,1,...,nrows-1</em> (if <em>indexFrom1</em>==<b>false</b>) or over <em>1,2,...,nrows</em> (if <em>indexFrom1</em>==<b>true</b>).
     *  \return  The permuted matrix <em>B</em>.
     *  \remark  This member function invokes permute_csc_rows() and then sort_csc_elements_in_each_column(),
     *           so that the return sparse matrix has elements in each column sorted with respect to row indices.
     *  \remark  This member function resembles <em>B</em>(<em>perm</em>,:)=<em>A</em> in <b>OCTAVE</b>/<b>MATLAB</b>.
     */
    // permute rows of A, such that
    //    B(perm[i-1],j)   = A(i,j) for i=1,...,nrows and j=1,...,ncols (if indexFrom1==true), or
    //    B(perm[i-1]+1,j) = A(i,j) for i=1,...,nrows and j=1,...,ncols (if indexFrom1==false), where
    //    A is *this and B is the return sparse matrix
    // the indices perm[0,...,nrows-1] range over 0,1,...,nrows-1 (if indexFrom1==false) or over 1,2,...,nrows (if indexFrom1==true)
    // this member function invokes permute_csc_rows() and then sort_csc_elements_in_each_column(),
    //    so that the return sparse matrix has elements in each column sorted with respect to row indices
    // this member function resembles B(perm,:) = A in OCTAVE/MATLAB
    SparseMatrix permuteRows(const mkIndex *perm, bool indexFrom1=true) const;

    //! Permute rows and columns of <em>*this</em>.
    /*! Let <em>A</em> be <em>*this</em> and <em>B</em> is the permuted matrix.
     *  <ul>
     *  <li>  If <em>indexFrom1</em>==<b>true</b>, <em>B</em>(<em>rperm</em>[<em>i-1</em>],<em>cperm</em>[<em>j-1</em>])               = <em>A</em>(<em>i,j</em>) for <em>i=1,...,nrows</em> and <em>j=1,...,ncols</em>.
     *  <li>  If <em>indexFrom1</em>==<b>false</b>, <em>B</em>(<em>rperm</em>[<em>i-1</em>]+<em>1,cperm</em>[<em>j-1</em>]+<em>1</em>) = <em>A</em>(<em>i,j</em>) for <em>i=1,...,nrows</em> and <em>j=1,...,ncols</em>.
     *  </ul>
     *  \param  rperm[] is the row permutation array.
     *  \param  cperm[] is the column permutation array.
     *  \param  indexFrom1 tells whether the indices start from <em>1</em> or <em>0</em>.
     *
     *  The indices <em>rperm</em>[<em>0,...,nrows-1</em>] range over <em>0,1,...,nrows-1</em> (if <em>indexFrom1</em>==<b>false</b>) or over <em>1,2,...,nrows</em> (if <em>indexFrom1</em>==<b>true</b>). \n
     *  The indices <em>cperm</em>[<em>0,...,ncols-1</em>] range over <em>0,1,...,ncols-1</em> (if <em>indexFrom1</em>==<b>false</b>) or over <em>1,2,...,ncols</em> (if <em>indexFrom1</em>==<b>true</b>).
     *  \remark  If <em>cperm</em>==<tt>NULL</tt>, then symmetric permutation will be applied (i.e., <em>rperm</em>[] will also be used as <em>cperm</em>[]).
     *  \remark  This member function resembles <em>B</em>(<em>rperm,cperm</em>)=<em>A</em> in <b>OCTAVE</b>/<b>MATLAB</b>.
     */
    // permute rows and columns of A, such that
    //    B(rperm[i-1],cperm[j-1])     = A(i,j) for i=1,...,nr and j=1,...,nc (if indexFrom1==true), or
    //    B(rperm[i-1]+1,cperm[j-1]+1) = A(i,j) for i=1,...,nr and j=1,...,nc (if indexFrom1==false), where
    //    A is the input matrix *this and B is the return output matrix
    // if cperm==NULL, then symmetric permutation will be applied (i.e., rperm[] will also be used for cperm[])
    // the indices rperm[0,...,nrows-1] range over 0,1,...,nrows-1 (if indexFrom1==false) or over 1,2,...,nrows (if indexFrom1==true)
    // the indices cperm[0,...,ncols-1] range over 0,1,...,ncols-1 (if indexFrom1==false) or over 1,2,...,ncols (if indexFrom1==true)
    // this member function resembles B(rperm,cperm) = A in OCTAVE/MATLAB
    SparseMatrix permuteRowsAndColumns(const mkIndex *rperm, const mkIndex *cperm=NULL, bool indexFrom1=true) const;

    #ifdef USE_MEX
    // convert *mxArray to SparseMatrix
    static SparseMatrix mxArray2SparseMatrix(const mxArray *mA);
    // convert SparseMatrix to *mxArray
    static mxArray *SparseMatrix2mxArray(const SparseMatrix &A);
    #endif
};


//! Return <em>f*A</em>.
/*! This operator (eventually) invokes vector_scalar_operation() to copy scaled nonzero elements.
 */
// sparse scalar - matrix multiplication
// return f*A
// this operator (eventually) invokes vector_scalar_operation() to copy scaled nonzero elements
SparseMatrix operator*(Real f, const SparseMatrix &A);


//! Display sparse matrix <em>A</em> with operator <<.
/*! In other words, print out sparse matrix <em>A</em> to <b>ostream</b> <em>s</em>.
 *  \remark  This routine invokes print_spmatrix().
 */
// display sparse matrix A with operator <<
// in other words, print out sparse matrix A to ostream s
// this routine invokes print_spmatrix()
std::ostream& operator<<(std::ostream &s, const SparseMatrix &A);


#ifdef USE_NAMESPACE
}  // end of namespace MATKIT
#endif

#endif  // end of #ifdef SPARSE_MATRIX_H
