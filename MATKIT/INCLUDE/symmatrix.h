#ifndef SYMMATRIX_H
#define SYMMATRIX_H
/*! \file symmatrix.h
 *  The class SymmetricMatrix defines a symmetric matrix represented by
 *  an array SymmetricMatrix::store[] ((with elements of the symmetric matrix in the lower triangular part stored columnwise in packed form)
 *  and the size SymmetricMatrix::nrows by SymmetricMatrix::nrows.
 *
 *  <ul>
 *  <li>  The <em>i</em>th element in the first column of the symmetric matrix
 *        is stored as SymmetricMatrix::store[<em>i-1</em>] for <em>i=1,...,</em>SymmetricMatrix::nrows.
 *  <li>  The <em>i</em>th element in the second column of the symmetric matrix
 *        is stored as SymmetricMatrix::store[(<em>i-2</em>)+ SymmetricMatrix::nrows] for <em>i=2,...,</em>SymmetricMatrix::nrows.
 *  <li>  In general, the (<em>i,j</em>) entry of the symmetric matrix is stored as
 *        SymmetricMatrix::store[(<em>i-j</em>) + (<em>j-1</em>)*(<em>2</em>* SymmetricMatrix::nrows-<em>j+2</em>)/<em>2</em>] for <em>1<=j<=i<=</em>SymmetricMatrix::nrows.
 *  </ul>
 */

#include <iostream>  // for ostream, etc.
#include "matkitdef.h"


#ifdef USE_NAMESPACE
namespace MATKIT {
#endif


class Vector;
class Matrix;
class SparseMatrix;

/*! \brief Class SymmetricMatrix defines a symmetric matrix and the arithmetic operations on it,
 *         and also providing some static functions for symmetric matrix generation.
 */
class SymmetricMatrix {
protected:
    //! The symmetric matrix <em>*this</em> has elements in the lower triangular part stored columnwise in packed form in <em>store</em>[].
    /*! To be precise, let <em>A</em> be the matrix <em>*this</em>. Then
     *  <em>A</em>(<em>i,j</em>) is stored in <em>store</em>[(<em>i-j</em>) + (<em>j-1</em>)*(<em>2*nrows-j+2</em>)/<em>2</em>] for <em>1<=j<=i<=nrows</em>.
     */
    // the symmetric matrix *this is stored columnwise in packed form in store[]
    // to be precise, let A be the matrix *this
    // then A(i,j) is stored in store[(i-j) + (j-1)*(2*nrows-j+2)/2] for 1<=j<=i<=nrows
    Real *store;
    //! Number of rows of <em>*this</em>, which is also the number of columns since the matrix is square.
    // number of rows of *this, which is also the number of columns since the matrix is square
    mkIndex nrows;

public:
    // some static functions
    //! Return an <em>n</em>-by-<em>n</em> (symmetric) matrix of zeros.
    // return an n-by-n (symmetric) matrix of zeros
    static SymmetricMatrix zeros(mkIndex n);
    //! Return an <em>n</em>-by-<em>n</em> (symmetric) matrix of ones.
    // return an n-by-n (symmetric) matrix of ones
    static SymmetricMatrix ones(mkIndex n);
    //! Return an <em>n</em>-by-<em>n</em> (symmetric) matrix with ones on the main diagonal and zeros elsewhere.
    // return an n-by-n (symmetric) matrix with ones on the main diagonal and zeros elsewhere
    static SymmetricMatrix eye(mkIndex n);
    //! Return an <em>n</em>-by-<em>n</em> symmetric matrix containing pseudo-random values drawn from a uniform distribution on the interval <em>[a,b]</em>.
    // return an n-by-n symmetric matrix containing pseudo-random values drawn from a uniform distribution on the interval [a,b]
    static SymmetricMatrix random(mkIndex n, Real a=0.0, Real b=1.0);


    // constructors
    //! A constructor for an empty (symmetric) matrix.
    // a constructor for an empty (symmetric) matrix
    SymmetricMatrix();
    //! A constructor for an <em>n</em>-by-<em>n</em> symmetric matrix of zeros.
    // a constructor for an n-by-n symmetric matrix of zeros
    SymmetricMatrix(mkIndex n);
    //! A constructor for an <em>n</em>-by-<em>n</em> symmetric matrix of with elements stored columnwise in packed form in <em>store0</em>[].
    /*! \remark  It performs a shallow copy of <em>store0</em>[] and the memory will be freed by the destructor (i.e. `<em>delete</em> [] <em>store0</em>').
     */
    // a constructor for an n-by-n symmetric matrix of with elements stored columnwise in packed form in store0[]
    // note that it performs a shallow copy of store0 and the memory will be freed by the destructor (i.e. `delete [] store0')
    SymmetricMatrix(mkIndex n, Real *store0);
    //! A copy constructor.
    /*! \remark  This constructor invokes memory_xcopy().
     */
    // a copy constructor
    // this constructor invokes memory_xcopy()
    SymmetricMatrix(const SymmetricMatrix &B);

    //! A destructor. Memory allocated, with the address pointed to by SymmetricMatrix::store, will be freed.
    // a destructor; memory allocated, with the address pointed to by SymmetricMatrix::store, will be freed
    virtual ~SymmetricMatrix();


    // get the information of *this
    //! The number of rows of <em>*this</em>.
    // the number of rows of *this
    mkIndex Nrows() const { return nrows; }
    //! The number of columns of <em>*this</em>.
    // the number of columns of *this
    mkIndex Ncols() const { return nrows; }
    //! A <b>const</b> pointer to the storage for elements of <em>*this</em>.
    // a const pointer to the storage for elements of *this
    const Real *Store() const { return store; }

    //! Convert <em>*this</em> to a general matrix.
    /*! This member function invokes symmetric_to_general().
     */
    // convert *this to a general matrix
    // this member function invokes symmetric_to_general()
    Matrix toGeneral() const;

    //! Resize <em>*this</em> to be of dimensions <em>nrc</em>-by-<em>nrc</em>.
    // resize *this to be of dimensions nrc-by-nrc
    void resize(mkIndex nrc);

    //! Set <em>*this</em> as an <em>nrc</em>-by-<em>nrc</em> symmetric matrix stored in <em>store0</em>[] columnwise in packed form.
    /*! It performs a shallow copy of <em>store0</em>[] and the memory will be freed by the destructor (i.e. `<em>delete</em> [] <em>store0</em>').
     */
    // set *this as an nrc-by-nrc symmetric matrix stored in store0[] columnwise in packed form
    // it performs a shallow copy of store0[] and the memory will be freed by the destructor (i.e. `delete [] store0')
    void set(mkIndex nrc, Real *store0);


    //! A vector formed by column <em>j</em> (or equivalently row <em>j</em>).
    /*! This member function invokes elements_of_symmatrix().
     */
    // a vector formed by column j (or equivalently row j)
    // this member function invokes elements_of_symmatrix()
    Vector column(mkIndex j) const;
    //! A vector formed by row <em>i</em> (or equivalently column <em>i</em>).
    /*! This member function invokes elements_of_symmatrix().
     */
    // a vector formed by row i (or equivalently column i)
    // this member function invokes elements_of_symmatrix()
    Vector row(mkIndex i) const;

    //! A submatrix formed by columns <em>j1,...,j2</em>.
    /*! This member function invokes submatrix_of_symmatrix().
     */
    // a submatrix formed by columns j1,...,j2
    // this member function invokes submatrix_of_symmatrix()
    Matrix columns(mkIndex j1, mkIndex j2) const;
    //! A submatrix formed by rows <em>i1,...,i2</em>.
    /*! This member function invokes submatrix_of_symmatrix().
     */
    // a submatrix formed by rows i1,...,i2
    // this member function invokes submatrix_of_symmatrix()
    Matrix rows(mkIndex i1, mkIndex i2) const;
    //! A submatrix formed by elements in rows <em>i1,...,i2</em> and columns <em>j1,...,j2</em>.
    /*! This member function invokes submatrix_of_symmatrix().
     */
    // a submatrix formed by elements in rows i1,...,i2 and columns j1,...,j2
    // this member function invokes submatrix_of_symmatrix()
    Matrix subMatrix(mkIndex i1, mkIndex i2, mkIndex j1, mkIndex j2) const;
    //! A submatrix as a symmetric matrix formed by elements in rows <em>k1,...,k2</em> and columns <em>k1,...,k2</em> of *this
    /*! This member function invokes submatrix_of_symmatrix().
     */
    // a submatrix as a symmetric matrix formed by elements in rows k1,...,k2 and columns k1,...,k2 of *this
    // this member function invokes submatrix_of_symmatrix()
    SymmetricMatrix subMatrix(mkIndex k1, mkIndex k2) const;


    //! A vector formed by the elements in the <em>k</em>th diagonal of <em>*this</em> in order.
    /*! This member function invokes symmatrix_diag().
     */
    // a vector formed by the elements in the k-th diagonal of *this in order
    // this member function invokes symmatrix_diag()
    Vector diag(mkSignedIndex k=0) const;

    //! The transpose of <em>*this</em>.
    // the transpose of *this
    SymmetricMatrix transpose() const { return *this; }


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
    bool operator==(const SymmetricMatrix &B) const {
        if (nrows != B.nrows)
            return false;
        return are_elements_equal_elements(nrows*(nrows+1u)/2u, Store(), B.Store());
    }
    //! Return <b>false</b> if <em>A==B</em> mathematically, and otherwise return <b>true</b>, where <em>A</em> is <em>*this</em>.
    /*! \remark  This operator (eventually) invokes are_elements_equal_elements().
     */
    // return false if A == B mathematically, and otherwise return true, where A is *this
    // this operator (eventually) invokes are_elements_equal_elements()
    bool operator!=(const SymmetricMatrix &B) const {
        return !(*this==B);
    }


    // operator =
    //! Return a <b>const</b> reference to <em>A</em> after setting <em>A</em>(<em>i,j</em>)=<em>f</em> for <em>i,j=1,...,nrows</em>, where <em>A</em> is <em>*this</em>.
    // return a const reference to A after setting A(i,j)=f for i,j=1,...,nrows, where A is *this
    const SymmetricMatrix& operator=(Real f);
    //! Return a <b>const</b> reference to <em>A</em> after setting <em>A=B</em>, where <em>A</em> is <em>*this</em>.
    /*! This operator invokes memory_xcopy().
     */
    // return a const reference to A after setting A=B, where A is *this
    // this operator invokes memory_xcopy()
    const SymmetricMatrix& operator=(const SymmetricMatrix &B);


    // operator -
    //! Return -<em>A</em>, where <em>A</em> is <em>*this</em>.
    /*! \remark  This operator invokes vector_scalar_operation().
     */
    // return -A, where A is *this
    // this operator invokes vector_scalar_operation()
    SymmetricMatrix operator-() const;


    // operators +,-,*,/,+=,-=,*-,/* a real number
    //! Return <em>B</em> with <em>B</em>(<em>i,j</em>)=<em>A</em>(<em>i,j</em>)+<em>f</em> for <em>i,j=1,...,nrows</em>, where <em>A</em> is <em>*this</em>.
    /*! \remark  This operator invokes vector_scalar_operation().
     */
    // return B with B(i,j)=A(i,j)+f for i,j=1,...,nrows, where A is *this
    // this operator invokes vector_scalar_operation()
    SymmetricMatrix operator+(Real f) const;
    //! Return <em>B</em> with <em>B</em>(<em>i,j</em>)=<em>A</em>(<em>i,j</em>)-<em>f</em> for <em>i,j=1,...,nrows</em>, where <em>A</em> is <em>*this</em>.
    /*! \remark  This operator invokes vector_scalar_operation().
     */
    // return B with B(i,j)=A(i,j)-f for i,j=1,...,nrows, where A is *this
    // this operator invokes vector_scalar_operation()
    SymmetricMatrix operator-(Real f) const;
    //! Return <em>B</em> with <em>B</em>(<em>i,j</em>)=<em>A</em>(<em>i,j</em>)*<em>f</em> for <em>i,j=1,...,nrows</em>, where <em>A</em> is <em>*this</em>.
    /*! \remark  This operator invokes vector_scalar_operation().
     */
    // return B with B(i,j)=A(i,j)*f for i,j=1,...,nrows, where A is *this
    // this operator invokes vector_scalar_operation()
    SymmetricMatrix operator*(Real f) const;
    //! Return <em>B</em> with <em>B</em>(<em>i,j</em>)=<em>A</em>(<em>i,j</em>)/<em>f</em> for <em>i,j=1,...,nrows</em>, where <em>A</em> is <em>*this</em>.
    /*! \remark  This operator invokes vector_scalar_operation().
     */
    // return B with B(i,j)=A(i,j)/f for i,j=1,...,nrows, where A is *this
    // this operator invokes vector_scalar_operation()
    SymmetricMatrix operator/(Real f) const;
    //! Return a <b>const</b> reference to <em>A</em> after computing <em>A</em>(<em>i,j</em>)+=<em>f</em> for <em>i,j=1,...,nrows</em>, where <em>A</em> is <em>*this</em>.
    /*! \remark  This operator invokes vector_scalar_operation().
     */
    // return a const reference to A after computing A(i,j)+=f for i=1,...,nrows and j=1,...,ncols, where A is *this
    // this operator invokes vector_scalar_operation()
    const SymmetricMatrix& operator+=(Real f);
    //! Return a <b>const</b> reference to <em>A</em> after computing <em>A</em>(<em>i,j</em>)-=<em>f</em> for <em>i,j=1,...,nrows</em>, where <em>A</em> is <em>*this</em>.
    /*! \remark  This operator (eventually) invokes vector_scalar_operation().
     */
    // return a const reference to A after computing A(i,j)-=f for i=1,...,nrows and j=1,...,ncols, where A is *this
    // this operator (eventually) invokes vector_scalar_operation()
    const SymmetricMatrix& operator-=(Real f);
    //! Return a <b>const</b> reference to <em>A</em> after computing <em>A</em>(<em>i,j</em>)*=<em>f</em> for <em>i,j=1,...,nrows</em>, where <em>A</em> is <em>*this</em>.
    /*! \remark  This operator invokes vector_scalar_operation().
     */
    // return a const reference to A after computing A(i,j)*=f for i=1,...,nrows and j=1,...,ncols, where A is *this
    // this operator invokes vector_scalar_operation()
    const SymmetricMatrix& operator*=(Real f);
    //! Return a <b>const</b> reference to <em>A</em> after computing <em>A</em>(<em>i,j</em>)/=<em>f</em> for <em>i,j=1,...,nrows</em>, where <em>A</em> is <em>*this</em>.
    /*! \remark  This operator invokes vector_scalar_operation().
     */
    // return a const reference to A after computing A(i,j)/=f for i=1,...,nrows and j=1,...,ncols, where A is *this
    // this operator invokes vector_scalar_operation()
    const SymmetricMatrix& operator/=(Real f);


    // symmetric matrix - symmetric matrix addition/subtraction
    //! Return <em>A+B</em>, where <em>A</em> is <em>*this</em>.
    /*! \remark  This operator invokes elementwise_addition().
     */
    // return A+B, where A is *this
    // this operator invokes elementwise_addition()
    SymmetricMatrix operator+(const SymmetricMatrix &B) const;
    //! Return <em>A-B</em>, where <em>A</em> is <em>*this</em>.
    /*! \remark  This operator invokes elementwise_addition().
     */
    // return A-B, where A is *this
    // this operator invokes elementwise_addition()
    SymmetricMatrix operator-(const SymmetricMatrix &B) const;
    //! Return a <b>const</b> reference to <em>A</em> after computing <em>A=A+B</em>, where <em>A</em> is <em>*this</em>.
    /*! \remark  This operator invokes elementwise_addition().
     */
    // return a const reference to A after computing A = A+B, where A is *this
    // this operator invokes elementwise_addition()
    const SymmetricMatrix& operator+=(const SymmetricMatrix &B);
    //! Return a <b>const</b> reference to <em>A</em> after computing <em>A=A-B</em>, where <em>A</em> is <em>*this</em>.
    /*! \remark  This operator invokes elementwise_addition().
     */
    // return a const reference to A after computing A = A-B, where A is *this
    // this operator invokes elementwise_addition()
    const SymmetricMatrix& operator-=(const SymmetricMatrix &B);

    // symmetric matrix - matrix addition/subtraction
    //! Return <em>A+B</em>, where <em>A</em> is <em>*this</em>.
    /*! \remark  This operator invokes matrix_symmatrix_addition().
     */
    // return A+B, where A is *this
    // this operator invokes matrix_symmatrix_addition()
    Matrix operator+(const Matrix &B) const;
    //! Return <em>A-B</em>, where <em>A</em> is <em>*this</em>.
    /*! \remark  This operator invokes matrix_symmatrix_addition().
     */
    // return A-B, where A is *this
    // this operator invokes matrix_symmatrix_addition()
    Matrix operator-(const Matrix &B) const;


    // symmetric matrix - vector multiplication
    //! Return <em>A*v</em>, where <em>A</em> is <em>*this</em>.
    /*! \remark  This operator invokes matrix_vector_multiplication().
     */
    // return A*v, where A is *this
    // this member function invokes symmatrix_vector_multiplication()
    Vector operator*(const Vector &v) const;

    // symmetric matrix - symmetric matrix multiplication
    //! Return <em>A*B</em>, where <em>A</em> is <em>*this</em>.
    /*! \remark  This operator invokes symmatrix_symmatrix_multiplication().
     */
    // return A*B, where A is *this
    // this member function invokes symmatrix_symmatrix_multiplication()
    Matrix operator*(const SymmetricMatrix &B) const;

    // symmetric matrix - matrix multiplication
    //! Return <em>A*B</em>, where <em>A</em> is <em>*this</em>.
    /*! \remark  This operator invokes symmatrix_matrix_multiplication().
     */
    // this member function invokes symmatrix_matrix_multiplication()
    Matrix operator*(const Matrix &B) const;

    //! Symmetrically permute rows and columns of <em>*this</em>.
    /*! Let <em>A</em> be <em>*this</em> and <em>B</em> is the permuted (symmetric) matrix.
     *  <ul>
     *  <li>  If <em>indexFrom1</em>==<b>true</b>, <em>B</em>(<em>perm</em>[<em>i-1</em>],<em>perm</em>[<em>j-1</em>]) = <em>A</em>(<em>i,j</em>) for <em>i,j=1,...,nrows</em>.
     *  <li>  If <em>indexFrom1</em>==<b>false</b>, <em>B</em>(<em>perm</em>[<em>i-1</em>]+<em>1,perm</em>[<em>j-1</em>]+<em>1</em>) = <em>A</em>(<em>i,j</em>) for <em>i,j=1,...,nrows</em>.
     *  </ul>
     *  \param  perm[]      is the permutation array.
     *  \param  indexFrom1  tells whether the indices start from <em>1</em> or <em>0</em>.
     *
     *  The indices <em>perm</em>[<em>0,...,nrows-1</em>] range over <em>0,1,...,nrows-1</em> (if <em>indexFrom1</em>==<b>false</b>) or over <em>1,2,...,nrows</em> (if <em>indexFrom1</em>==<b>true</b>).
     *
     *  \remark  This member function invokes permute_symmatrix_rows_and_columns().
     *  \remark  This member function resembles <em>B</em>(<em>perm,perm</em>)=<em>A</em> in <b>OCTAVE</b>/<b>MATLAB</b>.
     */
    // symmetrically permute rows and columns of A, such that
    //    B(perm[i-1]  ,perm[j-1]  ) = A(i,j) for i,j=1,...,nrows (if indexFrom1==true), or
    //    B(perm[i-1]+1,perm[j-1]+1) = A(i,j) for i,j=1,...,nrows (if indexFrom1==false),
    // where A is the input matrix *this and B is the returned output matrix
    // the indices perm[0,...,nrows-1] range over 0,1,...,nrows-1 (indexFrom1==false) or over 1,2,...,nrows (indexFrom1==true)
    // this member function resembles B(perm,perm) = A in OCTAVE/MATLAB
    // this member function invokes permute_symmatrix_rows_and_columns()
    SymmetricMatrix permuteRowsAndColumns(const mkIndex *perm, bool indexFrom1=true) const;
};


// symmetric scalar - matrix operations
//! Return <em>f*A</em>.
/*! This operator (eventually) invokes vector_scalar_operation().
 */
// return f*A
// this operator (eventually) invokes vector_scalar_operation()
SymmetricMatrix operator*(Real h, const SymmetricMatrix &A);

//! Return <em>B</em> with <em>B</em>(<em>i,j</em>)=<em>f+A</em>(<em>i,j</em>) for <em>i,j=1,...,A.nrows</em>.
/*! This operator (eventually) invokes vector_scalar_operation().
 */
// return B with B(i,j)=f+A(i,j) for i,j=1,...,A.nrows
// this operator (eventually) invokes vector_scalar_operation()
SymmetricMatrix operator+(Real h, const SymmetricMatrix &A);

//! Return <em>B</em> with <em>B</em>(<em>i,j</em>)=<em>f-A</em>(<em>i,j</em>) for <em>i,j=1,...,A.nrows</em>.
/*! This operator (eventually) invokes vector_scalar_operation().
 */
// return B with B(i,j)=f-A(i,j) for i,j=1,...,A.nrows
// this operator (eventually) invokes vector_scalar_operation()
SymmetricMatrix operator-(Real h, const SymmetricMatrix &A);


//! Display symmetric matrix <em>A</em> with operator <<.
/*! In other words, print out symmetric matrix <em>A</em> to <b>ostream</b> <em>s</em>.
 *  \remark  This routine invokes print_symmatrix().
 */
// display symmetric matrix A with operator <<
// in other words, print out symmetric matrix A to ostream s
// this routine invokes print_symmatrix()
std::ostream& operator<<(std::ostream &s, const SymmetricMatrix &A);


#ifdef USE_NAMESPACE
}  // end of namespace MATKIT
#endif

#endif  // end of #ifdef MATRIX_H
