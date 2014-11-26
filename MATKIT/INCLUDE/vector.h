#ifndef VECTOR_H
#define VECTOR_H
/*! \file vector.h
 *  The class Vector defines a vector represented by
 *  an array Vector::store[] and the length Vector::length.
 *  The <em>i</em>th element of the vector is stored as Vector::store[<em>i-1</em>] for <em>i=1,...,</em>Vector::length.
 */

#include <iostream>  // for ostream, etc.
#include "matkitdef.h"


#ifdef USE_NAMESPACE
namespace MATKIT {
#endif


class Matrix;
class SymmetricMatrix;
class SparseMatrix;


/*! \brief Class Vector defines a vector and the arithmetic operations on it,
 *         and also providing some static functions for matrix generation.
 */
class Vector {
    // friend classes and functions
    #ifdef USE_MWBLAS
    // MWBLAS routines dgemv (double precision) and sgemv (single precision) for matrix - vector multiplication demand non-const access to the matrix but indeed only const access is required
    friend class Matrix;
    // MWBLAS routines dspmv (double precision) and sspmv (single precision) for symmetric matrix - vector multiplication demand non-const access to the symmetric matrix but indeed only const access is required
    friend class SymmetricMatrix;
    #endif

protected:
    //! The vector <em>*this</em> is stored in <em>store</em>[].
    /*! To be precise, <em>v</em>(<em>i</em>) is stored in <em>store</em>[<em>i-1</em>] for <em>i=1,...,length</em>, where <em>v</em> is <em>*this</em>.
     */
    // the vector *this is stored in store[]
    // to be precise, v(i) is stored in store[i-1] for i=1,...,length, where v is *this
    Real *store;
    //! Length of <em>*this</em>.
    // length of *this
    mkIndex length;

public:
    //! Read a vector from a plain text file <em>fileName</em>[] and return the result.
    /*! \param  fileName[]  is the name of the plain text file.
     *  \remark  This member function invokes vector_read(); see the comments therein for the file format.
     */
    // read a vector from a plain text file fileName[] and return the result
    // this member function invokes vector_read(); see the comments therein for the file format
    static Vector read(const char fileName[]);

    //! Write a vector <em>v</em> to a plain text file <em>fileName</em>[], where <em>v</em> is <em>*this</em>.
    /*! \param  fileName[]  is the name of the plain text file.
     *  \remark  This member function invokes vector_write(); see the comments therein for the file format and the meaning of the returned value.
     */
    // write a vector v to a plain text file fileName[], where v is *this
    // this member function invokes vector_write(); see the comments therein for the file format and the meaning of the returned value
    int write(const char fileName[]) {
        return vector_write(length, store, fileName);
    }


    // some static functions
    //! Return a vector of zeros of length <em>n</em>.
    // return a vector of zeros of length n
    static Vector zeros(mkIndex n);
    //! Return a vector of ones of length <em>n</em>.
    // return a vector of ones of length n
    static Vector ones(mkIndex n);
    //! Return a vector of length <em>n</em> containing pseudo-random values drawn from a uniform distribution on the interval [<em>a,b</em>].
    // return a vector of length n containing pseudo-random values drawn from a uniform distribution on the interval [a,b]
    static Vector random(mkIndex n, Real a=0.0, Real b=1.0);
    //! Concatenate two vectors <em>v</em> and <em>w</em>.
    /*! \remark This static function invokes vector_concatenate().
     */
    // concatenate two vectors v and w
    // this static function invokes vector_concatenate()
    static Vector concatenate(const Vector &v, const Vector &w);
    //! Inner product <em>v'*w</em>.
    /*! \remark This static function invokes vector_inner_product().
     */
    // inner product v'*w
    // this static function invokes vector_inner_product()
    static Real innerProduct(const Vector &v, const Vector &w);
    //! (Scaled) outer product <em>alp*v*w'</em> as a general matrix.
    /*! \remark This static function invokes vector_outer_product(<b>mkIndex</b>, <b>mkIndex</b>, Real, const Real*, const Real*, Real*&, bool).
     */
    // (scaled) outer product alp*v*w' as a general matrix
    // this static function invokes vector_outer_product(mkIndex, mkIndex, Real, const Real*, const Real*, Real*&, bool)
    static Matrix outerProduct(const Vector &v, const Vector &w, Real alp=1.0);
    //! (Scaled) outer product <em>alp*v*v'</em> as a symmetric matrix.
    /*! \remark This static function invokes vector_outer_product(<b>mkIndex</b>, Real, const Real*, Real*&, bool).
     */
    // (scaled) outer product alp*v*v' as a symmetric matrix
    // this static function invokes vector_outer_product(mkIndex, Real, const Real*, Real*&, bool)
    static SymmetricMatrix outerProduct(const Vector &v, Real alp=1.0);
    //! (Scaled) symmetric outer product <em>alp</em>*(<em>v*w'+w*v'</em>) as a symmetric matrix.
    /*! \remark This static function invokes vector_symmetric_outer_product().
     */
    // (scaled) symmetric outer product alp*(v*w'+w*v') as a symmetric matrix
    // this static function invokes vector_symmetric_outer_product()
    static SymmetricMatrix symmetricOuterProduct(const Vector &v, const Vector &w, Real alp=1.0);


    // constructors
    //! A constructor for an empty vector.
    // a constructor for an empty vector
    Vector();
    //! A constructor for a vector of length n of zeros.
    // a constructor for a vector of length n of zeros

    Vector(mkIndex n);
    //! A constructor for a vector <em>v</em> of length <em>n</em> with elements stored in <em>store0</em>[] in order.
    /*! \remark It performs a shallow copy of <em>store0</em>[] and the memory will be freed by the destructor (i.e. `<em>delete</em> [] <em>store0</em>').
     */
    // a constructor for a vector v of length n with elements stored in store0[] in order
    // note that it performs a shallow copy of store0[] and the memory will be freed by the destructor (i.e. `delete [] store0')
    Vector(mkIndex n, Real *store0);

    //! A copy constructor.
    /*! \remark This constructor invokes memory_xcopy().
     */
    // a copy constructor
    // this constructor invokes memory_xcopy()
    Vector(const Vector &v);

    //! A destructor. Memory allocated, with the address pointed to by Vector::store, will be freed.
    // a destructor; memory allocated, with the address pointed to by Vector::store, will be freed
    virtual ~Vector();


    // get the information of *this
    //! The number of rows of <em>*this</em>.
    // the number of rows of *this
    mkIndex Length() const { return length; }
    //! The number of columns of <em>*this</em>.
    // the number of columns of *this
    const Real *Store() const { return store; }

    //! Resize <em>*this</em> to be of length <em>n</em>
    // resize *this to be of length n
    void resize(mkIndex n);
    //! Set <em>*this</em> be a vector of length <em>n</em> and stored in store[]
    /*! \remark It performs a shallow copy of <em>store0</em>[] and the memory will be freed by the destructor (i.e. `<em>delete</em> [] <em>store0</em>').
     */
    // set *this as a vector of length n stored in store0[]
    // it performs a shallow copy of store0[] and the memory will be freed by the destructor (i.e. `delete [] store0')
    void set(mkIndex n, Real *store0);

    //! Return a (square) matrix of dimension <em>len</em>+|<em>k</em>| with the <em>k</em>th diagonal formed by the elements of <em>*this</em> in order, and otherwise zero.
    /*! \remark This member function invokes vector_to_diag().
     *  \remark This member function resembles <b>diag</b>(<em>v</em>) in <b>OCTAVE</b>/<b>MATLAB</b>.
     */
    // return a (square) matrix of dimension len+|k| with the k-th diagonal formed by the elements of *this in order, and otherwise zero
    // this member function invokes vector_to_diag()
    // this member function resembles diag(v) in OCTAVE/MATLAB
    Matrix diag(mkSignedIndex k=0) const;
    //! Return a (square) sparse matrix of dimension <em>len</em>+|<em>k</em>| with the <em>k</em>th diagonal formed by the elements of <em>*this</em> in order, and otherwise zero.
    /*! \remark This member function invokes vector_to_spdiag().
     *  \remark This member function resembles <b>spdiag</b>(<em>v</em>) in <b>OCTAVE</b>/<b>MATLAB</b>.
     */
    // return a (square) sparse matrix of dimension len+|k| with the k-th diagonal formed by the elements of *this in order, and otherwise zero
    // this member function invokes vector_to_spdiag()
    // this member function resembles spdiag(v) in OCTAVE/MATLAB
    SparseMatrix spdiag(mkSignedIndex k=0) const;

    //! A subvector formed by elements <em>i,...,j</em> of <em>*this</em>.
    /*! \remark This member function invokes memory_xcopy().
     */
    // a subvector formed by elements i,...,j of *this
    // this member function invokes memory_xcopy()
    Vector subVector(mkIndex i, mkIndex j) const;


    // norms
    //! Frobenius norm of <em>*this</em>, the same as 2-norm since it is a vector.
    // Frobenius norm of *this, the same as 2-norm since it is a vector
    Real normFrobenius() const;
    //! 2-norm of <em>*this</em>.
    // 2-norm of *this
    Real norm2() const;
    //! 1-norm of <em>*this</em>.
    // 1-norm of *this
    Real norm1() const;
    //! \f$\infty\f$-norm of <em>*this</em>.
    // infinity-norm of *this
    Real normInfinity() const;
    //! Sum of squared elements of <em>*this</em>.
    // sum of squared elements of *this
    Real sumSquare() const;


    // access v(i)
    //! Return a reference to <em>v</em>(<em>i</em>), i.e. can be used to write <em>v</em>(<em>i</em>), where <em>v</em> is <em>*this</em>.
    // return a reference to v(i), i.e. can be used to write v(i), where v is *this
    virtual Real& operator()(mkIndex i);
    //! Read <em>v</em>(<em>i</em>) via a const reference or pointer to <em>v</em>, where <em>v</em> is <em>*this</em>.
    // read v(i) via a const reference or pointer to v, where v is *this
    Real operator()(mkIndex i) const;


    // operators ==, !=
    //! Return <b>true</b> if <em>w==v</em> mathematically, and otherwise return <b>false</b>, where <em>w</em> is <em>*this</em>.
    /*! \remark This operator invokes are_elements_equal_elements().
     */
    // return true if w == v mathematically, and otherwise return false, where w is *this
    // this operator invokes are_elements_equal_elements()
    bool operator==(const Vector &v) const {
        if (length != v.Length())
            return false;
        return are_elements_equal_elements(length, Store(), v.Store());
    }
    //! Return <b>false</b> if <em>w==v</em> mathematically, and otherwise return <b>true</b>, where <em>w</em> is <em>*this</em>.
    /*! \remark This operator (eventually) invokes are_elements_equal_elements().
     */
    // return false if w == v mathematically, and otherwise return true, where w is *this
    // this operator (eventually) invokes are_elements_equal_elements()
    bool operator!=(const Vector &v) const {
        return !(*this==v);
    }


    // operator =
    //! Return a const reference to <em>v</em> after setting <em>v</em>(<em>i</em>)=<em>f</em> for i=1,...,length, where <em>v</em> is <em>*this</em>.
    // return a const reference to v after setting v(i)=f for i=1,...,length, where v is *this
    const Vector& operator=(Real f);
    //! Return a <b>const</b> reference to <em>w</em> after setting <em>w=v</em>, where <em>w</em> is <em>*this</em>.
    /*! This operator invokes memory_xcopy().
     */
    // return a const reference to w after setting w=v, where w is *this
    // this operator invokes memory_xcopy()
    const Vector& operator=(const Vector &v);


    // operator -
    //! Return -<em>w</em>, where <em>w</em> is <em>*this</em>.
    /*! This operator invokes vector_scalar_operation().
     */
    // return -w, where w is *this
    // this operator invokes vector_scalar_operation()
    Vector operator-() const;

    // operator +,-,*,/,+=,-=,*=,/= a real number
    //! Return <em>v</em> with <em>v</em>(<em>i</em>)=<em>w</em>(<em>i</em>)+<em>f</em> for <em>i=1,...,length</em>, where <em>w</em> is <em>*this</em>.
    /*! \remark This operator invokes vector_scalar_operation().
     */
    // return v with v(i)=w(i)+f for i=1,...,length, where w is *this
    // this operator invokes vector_scalar_operation()
    Vector operator+(Real f) const;
    //! Return <em>v</em> with <em>v</em>(<em>i</em>)=<em>w</em>(<em>i</em>)-<em>f</em> for <em>i=1,...,length</em>, where <em>w</em> is <em>*this</em>.
    /*! \remark This operator invokes vector_scalar_operation().
     */
    // return v with v(i)=w(i)-f for i=1,...,length, where w is *this
    // this operator invokes vector_scalar_operation()
    Vector operator-(Real f) const;
    //! Return <em>v</em> with <em>v</em>(<em>i</em>)=<em>w</em>(<em>i</em>)*<em>f</em> for <em>i=1,...,length</em>, where <em>w</em> is <em>*this</em>.
    /*! \remark This operator invokes vector_scalar_operation().
     */
    // return v with v(i)=w(i)*f for i=1,...,length, where w is *this
    // this operator invokes vector_scalar_operation()
    Vector operator*(Real f) const;
    //! Return <em>v</em> with <em>v</em>(<em>i</em>)=<em>w</em>(<em>i</em>)/<em>f</em> for <em>i=1,...,length</em>, where <em>w</em> is <em>*this</em>.
    /*! \remark This operator invokes vector_scalar_operation().
     */
    // return v with v(i)=w(i)/f for i=1,...,length, where w is *this
    // this operator invokes vector_scalar_operation()
    Vector operator/(Real f) const;
    //! Return a <b>const</b> reference to <em>w</em> after computing <em>w</em>(<em>i</em>)+=<em>f</em> for <em>i=1,...,length</em>, where <em>w</em> is <em>*this</em>.
    /*! \remark This operator invokes vector_scalar_operation().
     */
    // return a const reference to w after computing w(i)+=f for i=1,...,length, where w is *this
    // this operator invokes vector_scalar_operation()
    const Vector& operator+=(Real f);
    //! Return a <b>const</b> reference to <em>w</em> after computing <em>w</em>(<em>i</em>)-=<em>f</em> for <em>i=1,...,length</em>, where <em>w</em> is <em>*this</em>.
    /*! \remark This operator invokes vector_scalar_operation().
     */
    // return a const reference to w after computing w(i)-=f for i=1,...,length, where w is *this
    // this operator invokes vector_scalar_operation()
    const Vector& operator-=(Real f);
    //! Return a <b>const</b> reference to <em>w</em> after computing <em>w</em>(<em>i</em>)*=<em>f</em> for <em>i=1,...,length</em>, where <em>w</em> is <em>*this</em>.
    /*! \remark This operator invokes vector_scalar_operation().
     */
    // return a const reference to w after computing w(i)*=f for i=1,...,length, where w is *this
    // this operator invokes vector_scalar_operation()
    const Vector& operator*=(Real f);
    //! Return a <b>const</b> reference to <em>w</em> after computing <em>w</em>(<em>i</em>)/=<em>f</em> for <em>i=1,...,length</em>, where <em>w</em> is <em>*this</em>.
    /*! \remark This operator invokes vector_scalar_operation().
     */
    // return a const reference to w after computing w(i)/=f for i=1,...,length, where w is *this
    // this operator invokes vector_scalar_operation()
    const Vector& operator/=(Real f);


    // vector - vector addition/subtraction
    //! Return <em>w+v</em>, where <em>w</em> is <em>*this</em>.
    /*! \remark This operator invokes elementwise_addition().
     */
    // return w + v, where w is *this
    // this operator invokes elementwise_addition()
    Vector operator+(const Vector &v) const;
    //! Return <em>w-v</em>, where <em>w</em> is <em>*this</em>.
    /*! \remark This operator invokes elementwise_addition().
     */
    // return w - v, where w is *this
    // this operator invokes elementwise_addition()
    Vector operator-(const Vector &v) const;
    //! Return a <b>const</b> reference to <em>w</em> after computing <em>w=w+v</em>, where <em>w</em> is <em>*this</em>.
    /*! \remark This operator invokes elementwise_addition().
     */
    // return a const reference to w after computing w = w+v, where w is *this
    // this operator invokes elementwise_addition()
    const Vector& operator+=(const Vector &v);
    //! Return a <b>const</b> reference to <em>w</em> after computing <em>w=w-v</em>, where <em>w</em> is <em>*this</em>.
    /*! \remark This operator invokes elementwise_addition().
     */
    // return a const reference to w after computing w = w-v, where w is *this
    // this operator invokes elementwise_addition()
    const Vector& operator-=(const Vector &v);


    // vector - matrix multiplication
    //! Compute <em>v'=w'*A</em> and return the result <em>v</em>, where <em>w</em> is <em>*this</em>.
    /*! This operator invokes vector_matrix_multiplication().
     */
    // compute v' = w'*A and return the result v, where w is *this
    // this operator invokes vector_matrix_multiplication()
    Vector operator*(const Matrix &A) const;
    //! Compute <em>w'=w'*A</em> and return a <b>const</b> reference to <em>w</em>, where <em>w</em> is <em>*this</em>.
    /*! This operator (eventually) invokes vector_matrix_multiplication().
     */
    // compute w' = w'*A and return the const reference to w, where w is *this
    // this operator (eventually) invokes vector_matrix_multiplication()
    const Vector& operator*=(const Matrix &A);


    // vector - symmetric matrix multiplication
    //! Compute <em>v'=w'*A</em> and return the result <em>v</em>, where <em>w</em> is <em>*this</em>.
    /*! This operator (eventually) invokes symmatrix_vector_multiplication().
     */
    // compute v' = w'*A and return the result v, where w is *this
    // this operator (eventually) invokes symmatrix_vector_multiplication()
    Vector operator*(const SymmetricMatrix &A) const;
    //! Compute <em>w'=w'*A</em> and return a <b>const</b> reference to <em>w</em>, where <em>w</em> is <em>*this</em>.
    /*! This operator (eventually) invokes symmatrix_vector_multiplication().
     */
    // compute w' = w'*A and return the const reference to w, where w is *this
    // this operator (eventually) invokes symmatrix_vector_multiplication()
    const Vector& operator*=(const SymmetricMatrix &A);


    // vector - sparse matrix multiplication
    //! Compute <em>v'=w'*A</em> and return the result <em>v</em>, where <em>w</em> is <em>*this</em>.
    /*! This operator invokes vector_spmatrix_multiplication().
     */
    // compute v' = w'*A and return the result v, where w is *this
    // this operator invokes vector_spmatrix_multiplication()
    Vector operator*(const SparseMatrix &A) const;
    //! Compute <em>w'=w'*A</em> and return a <b>const</b> reference to <em>w</em>, where <em>w</em> is <em>*this</em>.
    /*! This operator (eventually) invokes vector_matrix_multiplication().
     */
    // compute w' = w'*A and return the const reference to w, where w is *this
    // this operator (eventually) invokes vector_spmatrix_multiplication()
    const Vector& operator*=(const SparseMatrix &A);


    //! Return a vector with elements as those of <em>*this</em> reversed
    /*! To be precise, return <em>v</em> such that <em>v</em>(<em>i</em>)=<em>w</em>(<em>length-i+1</em>) for <em>i=1,...,length</em>, where <em>w</em> is <em>*this</em>.
     */
    // return a vector with elements as those of *this reversed
    // to be precise, return v such that v(i) = w(length-i+1) for i=1,...,length, where w is *this
    Vector reverse() const;


    //! Modified Gram-Schmidt, assuming input <em>Q</em> consists of orthonormal columns.
    /*! On exit, <em>*this</em> is orthogonal to columns of <em>Q</em>. \n
     *  If <em>normalize</em>==<b>true</b>, then on exist, <em>*this</em> will be normalized to have unit length, i.e. <em>this->Norm2</em>()==<em>1.0</em>.
     */
    // modified Gram-Schmidt, assuming input Q consists of orthonormal columns
    // on exit, *this is orthogonal to columns of Q
    // if normalize==true, then on exist, *this will be normalized to have unit length, i.e. this->Norm2()==1.0
    void modifiedGramSchmidt(const Matrix &Q, bool normalize=false);
    //! Modified Gram-Schmidt to orthogonalize against selected columns of <em>Q</em>, indicated by the bool array <em>flag</em>[] of length <em>Q.Ncols</em>().
    /*! On exit, <em>*this</em> is orthogonal to <em>Q</em>(:,<em>j</em>) for <em>flag</em>[<em>j-1</em>]==<b>true</b>. \n
     *  It is assumed that <em>Q</em>(:,<em>i</em>) and <em>Q</em>(:,<em>j</em>) are orthogonal to each other for <em>i</em>\f$\neq\f$<em>j</em>. \n
     *  If <em>normalize</em>==<b>true</b>, then on exist, <em>*this</em> will be normalized to have unit length, i.e. <em>this->Norm2</em>()==<em>1.0</em>.
     */
    // modified Gram-Schmidt to orthogonalize against selected columns of Q, indicated by the bool array flag[] of length Q.Ncols()
    // on exit, *this is orthogonal to Q(:,j) for flag[j-1]==true
    // it is assumed that Q(:,i) and Q(:,j) are orthogonal to each other for i != j
    // if normalize==true, then on exist, *this will be normalized to have unit length, i.e. this->Norm2()==1.0
    void modifiedGramSchmidt(const Matrix &Q, const bool *flag, bool normalize=false);

    #ifdef USE_MEX
    // convert *mxArray to Vector
    static Vector mxArray2Vector(const mxArray *mA);
    // convert Vector to *mxArray
    static mxArray *Vector2mxArray(const Vector &A);
    #endif
};


// scalar - vector operations
//! Return <em>w</em> with <em>w</em>(<em>i</em>)=<em>f+v</em>(<em>i</em>) for <em>i=1,...,v.length</em>.
/*! This operator (eventually) invokes vector_scalar_operation().
 */
// return w with w(i)=f+v(i) for i=1,...,v.length
// this operator (eventually) invokes vector_scalar_operation()
Vector operator+(Real f, const Vector &v);

//! Return <em>w</em> with <em>w</em>(<em>i</em>)=<em>f-v</em>(<em>i</em>) for <em>i=1,...,v.length</em>.
/*! This operator (eventually) invokes vector_scalar_operation().
 */
// return w with w(i)=f-v(i) for i=1,...,v.length
// this operator (eventually) invokes vector_scalar_operation()
Vector operator-(Real f, const Vector &v);

//! Return <em>f*v</em>.
/*! This operator (eventually) invokes vector_scalar_operation().
 */
// return f*v
// this operator (eventually) invokes vector_scalar_operation()
Vector operator*(Real f, const Vector &v);


//! Display vector <em>v</em> with operator <<.
/*! In other words, print out vector <em>v</em> to <b>ostream</b> <em>s</em>.
 *  \remark This routine invokes print_vector().
 */
// display vector v with operator <<
// in other words, print out vector v to ostream s
// this routine invokes print_vector()
std::ostream& operator<<(std::ostream &s, const Vector &v);


#ifdef USE_NAMESPACE
}  // end of namespace MATKIT
#endif

#endif  // end of #ifdef VECTOR_H
