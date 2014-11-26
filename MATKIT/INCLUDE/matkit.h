/*! \mainpage MATKIT
 *
 *
 *  \section sec_what What is MATKIT?
 *
 *  The <b>MATKIT</b> library is a collection of functions and classes written in <b>C++</b>
 *  for basic matrix computations.
 *
 *  The key features include:
 *  <ul>
 *  <li>  Several functions can optionally invoke <b>BLAS</b> routines to improve efficiency. \n
 *        (ps. This is particularly useful for multi-core parallelism if <b>Intel MKL BLAS</b> is readily available.)
 *  <li>  It supports both dense and sparse matrix computations (e.g. matrix - vector multiplication, with the matrix dense or sparse).
 *  <li>  It provides a set of functions and four basic classes Vector, Matrix, SymmetricMatrix, and SparseMatrix.
 *  <li>  Some functions resemble basic <b>OCTAVE</b>/<b>MATLAB</b> commands (e.g. <b>speye</b>(<em>n</em>), <b>full</b>(<em>A</em>), and <b>diag</b>(<em>A,k</em>)).
 *  <li>  A <b>MEX</b> interface is also provided for <b>OCTAVE</b>/<b>MATLAB</b> users.
 *  </ul>
 *
 *
 *  \section sec_data How are the matrices and vectors stored?
 *
 *  \subsection single_or_double Single precision or double precision?
 *
 *  The fundamental elements (i.e. entries in a vector or a matrix) are real numbers,
 *  which are represented as single precision (<b>float</b>) or
 *  double precision (<b>double</b>) floating point numbers.
 *
 *  By default, <b>double</b> is in use and defines the type <b>Real</b>.
 *  Alternatively, one may set <tt>USE_SINGLE</tt> for <b>float</b> by either `<tt>#define USE_SINGLE</tt>' in matkitdef.h or
 *  the compiler flag `<tt>-DUSE_SINGLE</tt>'.
 *  The latter way is used in the provided makefiles.
 *
 *  \subsection columnwise All are stored columnwise.
 *
 *  A vector <em>v</em> is represented by an array <em>store</em>[] and its length <em>length</em>.
 *      <ul>
 *      <li>  The <em>i</em>th element of <em>v</em>, also referred to as <em>v</em>(<em>i</em>), is stored in <em>store</em>[<em>i-1</em>] for <em>i=1,...,length</em>.
 *      <li>  Here <em>store</em> and <em>length</em> correspond to Vector::store and Vector::length of class Vector.
 *      </ul>
 *
 *  A general matrix <em>A</em> is represented by an array <em>store</em>[] (with elements of <em>A</em> stored columnwise) and its size <em>nrows</em>-by-<em>ncols</em>.
 *      <ul>
 *      <li>  The (<em>i,j</em>) entry of <em>A</em>, also referred to as <em>A</em>(<em>i,j</em>), is stored in <em>store</em>[(<em>i-1</em>)+(<em>j-1</em>)*<em>nrows</em>]
 *            for <em>i=1,...,nrows</em> and <em>j=1,...,ncols</em>.
 *      <li>  Here <em>store</em>, <em>nrows</em>, and <em>ncols</em> correspond to Matrix::store, Matrix::nrows, and Matrix::ncols of class Matrix.
 *      </ul>
 *
 *  A symmetric matrix <em>B</em> is represented by an array <em>store</em>[] (with elements of <em>B</em> in the lower triangular part stored columnwise in packed form) and its size <em>nrows</em>-by-<em>nrows</em>.
 *      <ul>
 *      <li>  The <em>i</em>th element in the first column of <em>B</em>, <em>B</em>(<em>i,1</em>), is stored in <em>store</em>[<em>i-1</em>] for <em>i=1,...,nrows</em>.
 *      <li>  The <em>i</em>th element in the second column of <em>B</em>, <em>B</em>(<em>i,2</em>), is stored in <em>store</em>[(<em>i-2</em>)+<em>nrows</em>] for <em>i=2,...,nrows</em>.
 *      <li>  In general, the (<em>i,j</em>) entry of <em>B</em>, <em>B</em>(<em>i,j</em>), is stored in
 *            <em>store</em>[(<em>i-j</em>) + (<em>j-1</em>)*(<em>2*nrows-j+2</em>)/<em>2</em>] for <em>1<=j<=i<=nrows</em>.
 *      <li>  Here <em>store</em> and <em>nrows</em> correspond to SymmetricMatrix::store and SymmetricMatrix::nrows of class SymmetricMatrix.
 *      </ul>
 *
 *  A sparse matrix <em>S</em> is stored in compressed sparse column (CSC) format in <em>store</em>[], <em>rowIndex</em>[], and <em>colPointer</em>[].
 *  In what follows we assume <em>colPointer</em>[<em>0</em>]==<em>0</em> for simplicity.
 *      <ul>
 *      <li>  The number of nonzero elements in <em>j</em>th column is <em>colPointer</em>[<em>j</em>]-<em>colPointer</em>[<em>j-1</em>].
 *      <li>  The nonzero elements in <em>j</em>th column are stored in <em>store</em>[<em>colPointer</em>[<em>j-1</em>],...,<em>colPointer</em>[<em>j</em>]-1].
 *      <li>  The corresponding row indices are stored in <em>rowIndex</em>[<em>colPointer</em>[<em>j-1</em>],...,<em>colPointer</em>[<em>j</em>]-1]. \n
 *            (ps. The index of <em>i</em>th row is <em>i-1</em>.)
 *      <li>  Here <em>store</em>, <em>rowIndex</em>, and <em>colPointer</em> correspond to SparseMatrix::store, SparseMatrix::rowIndex, SparseMatrix::colPointer of class SparseMatrix.
 *      </ul>
 *
 *  \remark The data structures in <b>OCTAVE</b>/<b>MATLAB</b> match those in <b>MATKIT</b>.
 *          So it is easy to import matrices, both dense and sparse, from <b>OCTAVE</b>/<b>MATLAB</b> to <b>MATKIT</b>
 *          and vice versa, via <b>MEX</b> files.
 *
 *  \subsection wo_class MATKIT without classes
 *
 *  There are a number of <b>MATKIT</b> functions without classes.
 *  Indeed, these functions use only intrinsic data types of <b>C++</b>, except for a few templates and
 *  defined <b>mkIndex</b> for indices and <b>Real</b> for real numbers.
 *  For example, as described above, a matrix is specified by two <b>mkIndex</b> numbers for number of rows and number of columns, and
 *  a pointer to the storage of <b>Real</b> numbers for storing elements of <em>A</em>.
 *
 *  We illustrate the function matrix_vector_multiplication() without classes as follows.
 *  This function computes the multiplication of a general matrix <em>A</em> and a vector <em>v</em>,
 *  and the output is a vector <em>w</em>.
 *  The interface is:
 *
 *  <tt>void matrix_vector_multiplication(mkIndex nrows, mkIndex ncols, const Real *sa, const Real *sv, Real *&sw, bool init0=true);</tt>
 *
 *  It takes four parameters for the input matrix and vector:
 *  <ul>
 *  <li>  <b>mkIndex</b> <em>nrows</em> is the number of rows.
 *  <li>  <b>mkIndex</b> <em>ncols</em> is the number of columns.
 *  <li>  <b>const</b> <b>Real</b> <em>sa</em>[] stores the elements of <em>A</em>. The value of <em>A</em>(<em>i,j</em>) is <em>sa</em>[(i-1)+(j-1)*nr] for <em>i=1,...,nrows</em> and <em>j=1,...,ncols</em>.
 *  <li>  <b>const</b> <b>Real</b> <em>sv</em>[] stores the elements of <em>v</em>. The value of <em>v</em>(<em>j</em>) is <em>sv</em>[<em>j-1</em>] for <em>j=1,...,ncols</em>.
 *  </ul>
 *  The output vector <em>w</em> is stored in <em>sw</em>[].
 *  The value of <em>w</em>(<em>i</em>) is <em>sw</em>[<em>i-1</em>] for <em>i=1,...,nrows</em>.
 *
 *  There is a <b>bool</b> flag <em>init0</em> whose default value is <b>true</b>.
 *  <ul>
 *  <li>  If <em>init0</em>==<b>true</b> or <em>sw</em>==<tt>NULL</tt>, then <em>sw</em>[<em>0,...,nrows-1</em>] are initialized as zeros and this function computes <em>w=A*v</em>. \n
 *  <li>  If <em>init0</em>==<b>false</b> and <em>sw</em>\f$\neq\f$<tt>NULL</tt>, then this function computes <em>w+=A*v</em>.
 *  </ul>
 *
 *  \remark matrix_vector_multiplication() invokes the level 2 <b>BLAS</b> routine <b>xGEMV</b> if <tt>USE_BLAS</tt> is defined.
 *
 *
 *  \section memory_management Memory management
 *
 *  <b>MATKIT</b> uses the <b>C++</b> operators <b>new</b> and <b>delete</b> to allocate and free memory.
 *  If you would like to use another memory management scheme (e.g. the <b>C</b> functions <b>malloc</b>() and <b>free</b>()),
 *  you can override the two <b>C++</b> operators.
 *
 *  For the <b>MATKIT</b> functions without classes, memory is allocated only for output storage and work space (i.e. temporary memory).
 *  In both cases, a user can allocate the memory and pass the pointers to the functions, in which case
 *  memory will not be allocated or freed in these functions.
 *  In other words, memory management can be handled completely outside these functions.
 *  On the other hand, for each output storage, a reference to a pointer is passed, and for each work space, a pointer is passed.
 *  In both cases, if the pointer points to <tt>NULL</tt>, then the required memory will be allocated inside the function, and
 *  the users do not have to worry about how much memory should be allocated.
 *  If work space is allocated inside a function, it will be freed on return.
 *  If output storage is allocated inside a function, then its address is stored in the pointer on return,
 *  since we passed a reference to the pointer.
 *
 *  In the above example matrix_vector_multiplication(), the parameter <em>sw</em> is passed as a reference to a pointer
 *  for the output vector <em>w</em>.
 *  If <em>sw</em>==<tt>NULL</tt> is given as input, then the memory of size <em>nrows</em> is allocated and pointed to by <em>sw</em>
 *  for the vector <em>w</em> which is of length <em>nrows</em>, and
 *  <em>sw</em> as output still points to the location of allocated memory.
 *
 *
 *  \section exception_handle Exception handling
 *
 *  The <b>MATKIT</b> functions without classes generally assume that the inputs are valid.
 *  For example, in the above function matrix_vector_multiplication(), it is assumed that
 *  the input pointer <em>sa</em> points to an address of memory allowing at least <em>nrows*ncols</em> <b>Real</b> numbers.
 *  On the other hand, some functions return an <b>int</b> flag which reflects the status.
 *  For example, if matrix_market_spmatrix_read() encounters a file read error, a negative number will be returned.
 *
 *  For class member functions, an error message is printed to Basics::err and
 *  the function Basics::quit() is called whenever an exception is catched.
 *  By default, Basics::err points to the <b>ostream</b> <b>std::err</b> and Basics::quit() calls the <b>C</b> function <b>exit</b>().
 *
 *
 *  \section mex_support MEX support for OCTAVE/MATLAB users
 *
 *  <b>MATKIT</b> can be used for <b>MEX</b> code with the compiler flag `<tt>MEX_FLAG=</tt>' (i.e. empty string) for
 *  <b>OCTAVE</b> or the flag `<tt>MEX_FLAG=USE_MATLAB</tt>' for <b>MATLAB</b>.
 *  The flag <tt>MEX_FLAG</tt> is included in the file `<tt>Make.inc</tt>'.
 *
 */
#ifndef MATKIT_H
#define MATKIT_H
/*! \file matkit.h
 *  This file is to include all <b>MATKIT</b> header files.
 */

#include "matkitdef.h"
#include "matkitfunc.h"

#include "basics.h"
#include "vector.h"
#include "matrix.h"
#include "symmatrix.h"
#include "spmatrix.h"

#endif  // end of #ifndef MATKIT_H
