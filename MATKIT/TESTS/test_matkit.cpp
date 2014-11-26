// This driver illustrates how to use MATKIT classes Vector, Matrix,
// SymmetricMatrix, and SparseMatrix.
// See also the OCTAVE/MATLAB script "test_matkit.m" for the corresponding
// operations tested.

#include <iostream>  // for cout, cerr, endl, etc. (under namespace std)
#include <time.h>    // for time_t, time, clock_t, clock, CLOCKS_PER_SEC, etc.
#include "matkit.h"

using std::cout;
using std::endl;

#ifdef USE_NAMESPACE
using namespace MATKIT;
#endif


int main(int argc, char *argv[]) {
    // this permutation arrays will be used for tests on permuting elements of a vector and rows/columns of a matrix
    const mkIndex perm2[] = { 2, 1 };        // indices starting from 1
    const mkIndex perm3[] = { 2, 3, 1 };     // indices starting from 1
    const mkIndex perm4[] = { 3, 2, 4, 1 };  // indices starting from 1
    // const mkIndex perm2x[] = { 1, 0 };        // indices starting from 0
    // const mkIndex perm3x[] = { 1, 2, 0 };     // indices starting from 0
    // const mkIndex perm4x[] = { 2, 1, 3, 0 };  // indices starting from 0


    ///////////////////////////////////////////////////////////////////////////
    //    machine precision
    ///////////////////////////////////////////////////////////////////////////
    cout << "### machine precision" << endl << endl;

    cout << "infinity  = " << Basics::infinity << endl;
    cout << "-infinity = " << -Basics::infinity << endl;
    cout << "eps = " << Basics::machineEpsilon << endl;
    cout << "1 + eps   - 1 = " << 1 + Basics::machineEpsilon - 1 << endl;
    cout << "1 + eps/2 - 1 = " << 1 + Basics::machineEpsilon/2 - 1 << endl;

    cout << endl << endl;


    ///////////////////////////////////////////////////////////////////////////
    //    vector routines
    ///////////////////////////////////////////////////////////////////////////
    cout << "### vector routines" << endl << endl;

    // vector construction
    Vector v(3);
    v(1) = 1.0;    v(2) = 2.0;    v(3) = 3.0;
    cout << "v =" << endl << v << endl;

    // subvector
    cout << "v.subVector(3,1) =" << endl << v.subVector(3,1) << endl;
    cout << "v.subVector(2,3) =" << endl << v.subVector(2,3) << endl;

    // vector norm
    cout << "||v||_2^2 = " << Basics::square(v.norm2()) << endl << endl;

    // for diag(v) as a general matrix
    cout << "v.diag() =" << endl << v.diag() << endl;
    // for diag(v,1) as a general matrix
    cout << "v.diag(1) =" << endl << v.diag(1) << endl;
    // for spdiag(v,-2) as a sparse matrix
    cout << "v.spdiag(-2) =" << endl << v.spdiag(-2) << endl;

    // vector - scalar operations
    cout << "-v=" << endl << -v << endl;
    cout << "3*v =" << endl << 3.0*v << endl;
    cout << "v/2 =" << endl << v/2.0 << endl;
    cout << "1+v =" << endl << 1.0+v << endl;
    cout << "v-2 =" << endl << v-2.0 << endl;
    cout << "v*=3;" << endl;
    cout << "v =" << endl << (v*=3.0) << endl;
    cout << "v-=2;" << endl;
    cout << "v =" << endl << (v-=2.0) << endl;
    cout << "v+=5;" << endl;
    cout << "v =" << endl << (v+=5.0) << endl;
    cout << "v/=3;" << endl;
    cout << "v =" << endl << (v/=3.0) << endl;

    // vector - vector addition/subtraction
    Vector w(3);
    w(1) = 1.0;    w(2) = -2.0;    w(3) = 4.0;
    cout << "w =" << endl << w << endl;
    cout << "v+w =" << endl << v+w << endl;
    cout << "v-w =" << endl << v-w << endl;
    cout << "v-=w;" << endl;
    cout << "v =" << endl << (v-=w) << endl;
    cout << "v+=w;" << endl;
    cout << "v =" << endl << (v+=w) << endl;

    // vector inner/outer products
    cout << "v'*w = " << Vector::innerProduct(v,w) << endl;
    // for v*v' as a symmetric matrix
    cout << "v*v' =" << endl << Vector::outerProduct(v) << endl;
    // for v*w' as a general matrix
    cout << "v*w' =" << endl << Vector::outerProduct(v,w) << endl;
    // for v*w'+w*v' as a symmetric matrix
    cout << "v*w'+w*v' = " << endl << Vector::symmetricOuterProduct(v,w) << endl;

    // vector == vector
    cout << "v == w?    " << (v==w) << endl;
    cout << "v != w?    " << (v!=w) << endl;
    cout << "v == v?    " << (v==v) << endl;
    cout << "v != v?    " << (v!=v) << endl;

    cout << endl << endl;


    ///////////////////////////////////////////////////////////////////////////
    //    matrix routines
    ///////////////////////////////////////////////////////////////////////////
    cout << "### matrix routines" << endl << endl;

    // function eye
    // for eye(3) as a general matrix
    cout << "Matrix::eye(2) =" << endl << Matrix::eye(2) << endl;
    // for eye(3,2) as a general matrix
    cout << "Matrix::eye(3,2) =" << endl << Matrix::eye(3,2) << endl;

    // function zeros
    // for zeros(2) as a general matrix
    cout << "Matrix::zeros(2) =" << endl << Matrix::zeros(2) << endl;
    // for zeros(3,2) as a general matrix
    cout << "Matrix::zeros(3,2) =" << endl << Matrix::zeros(3,2) << endl;

    // function ones
    // for ones(2) as a general matrix
    cout << "Matrix::ones(3) =" << endl << Matrix::ones(3) << endl;
    // for ones(2,3) as a general matrix
    cout << "Matrix::ones(3,1) =" << endl << Matrix::ones(3,1) << endl;

    // for sparse(ones(2,3)-eye(2,3))
    cout << "(Matrix::ones(2,3)-Matrix::eye(2,3)).toSparse() =" << endl << (Matrix::ones(2,3)-Matrix::eye(2,3)).toSparse() << endl;

    // matrix assignment
    // for A34 = [ 1  2  3  4; 5  6  7  8; 9  10  11  12 ]
    Matrix A34(3,4);
    Real val = 0.0;
    for (mkIndex i=1; i<=3; i++) {
        for (mkIndex j=1; j<=4; j++)
            A34(i,j) = (val+=1.0);
    }
    cout << "A34 =" << endl << A34 << endl;

    // submatrix of a general matrix
    mkIndex indices[] = { 3, 1 };
    // for A34(:,4) as a vector
    cout << "A34.column(4) =" << endl << A34.column(4) << endl;
    // for A34(2,:) as a vector
    cout << "A34.row(2) =" << endl << A34.row(2) << endl;
    // for A34(4:-1:2,:) as a matrix
    cout << "A34.columns(4,2) =" << endl << A34.columns(4,2) << endl;
    // for A34(:,1:3) as a matrix
    cout << "A34.columns(1,3) =" << endl << A34.columns(1,3) << endl;
    // for A34(:,[3,1]) as a matrix
    cout << "A34.columns([3,1]) =" << endl << A34.columns(2,indices) << endl;
    // for A34(2:3,:) as a matrix
    cout << "A34.rows(2,3) =" << endl << A34.rows(2,3) << endl;
    // for A34(3:-1:1,:) as a matrix
    cout << "A34.rows(3,1) =" << endl << A34.rows(3,1) << endl;
    // for A34([3,1],:) as a matrix
    cout << "A34.rows(2,[3,1]) =" << endl << A34.rows(2,indices) << endl;
    // for A34(3:-1:2,3:-1:1) as a matrix
    cout << "A34.subMatrix(3,2,3,1) =" << endl << A34.subMatrix(3,2,3,1) << endl;
    // for A34(1:2,2:4) as a matrix
    cout << "A34.subMatrix(1,2,2,4) =" << endl << A34.subMatrix(1,2,2,4) << endl;

    // matrix diagonal
    // for diag(A34) as a vector
    cout << "A34.diag()" << endl << A34.diag() << endl;
    // for diag(A34,2) as a vector
    cout << "A34.diag(2)" << endl << A34.diag(2) << endl;
    // for diag(A34,-1) as a vector
    cout << "A34.diag(-1)" << endl << A34.diag(-1) << endl;

    // matrix transpose
    // for A34'
    cout << "A34.transpose() =" << endl << A34.transpose() << endl;

    // norms of a general matrix
    cout << "||A34||_F^2 = " << Basics::square(A34.normFrobenius()) << endl;
    cout << "||A34||_1   = " << A34.norm1() << endl;
    cout << "||A34||_inf = " << A34.normInfinity() << endl << endl;

    // matrix xsum (and also rows and transpose)
    cout << "A24 = A34.rows(1,2);" << endl;
    Matrix A24 = A34.rows(1,2);
    cout << "A24 =" << endl << A24 << endl;
    cout << "A42 = A24.transpose();" << endl;
    Matrix A42 = A24.transpose();
    cout << "A42 =" << endl << A42 << endl;
    cout << "xsum(A42,A24) =" << endl << Matrix::xsum(A42,A24) << endl;

    // permutation of rows / columns of a general matrix
    cout << "perm3 = [ 2, 3, 1 ]" << endl;
    cout << "perm4 = [ 3, 2, 4, 1 ]" << endl;
    // cout << "perm3x = [ 1, 2, 0 ]" << endl;
    // cout << "perm4x = [ 2, 1, 3, 0 ]" << endl;
    cout << endl;
    cout << "A34.permuteColumns(perm4) =" << endl << A34.permuteColumns(perm4) << endl;
    // cout << "A34.permuteColumns(perm4x, false) =" << endl << A34.permuteColumns(perm4x, false) << endl;  // same result as above
    cout << "A34.permuteRows(perm3) =" << endl << A34.permuteRows(perm3) << endl;
    // cout << "A34.permuteRows(perm3x, false) =" << endl << A34.permuteRows(perm3x, false) << endl;  // same result as above
    cout << "A34.permuteRowsAndColumns(perm3, perm4) =" << endl << A34.permuteRowsAndColumns(perm3, perm4) << endl;
    // cout << "A34.permuteRowsAndColumns(perm3x, perm4x, false) =" << endl << A34.permuteRowsAndColumns(perm3x, perm4x, false) << endl;  // same result as above

    // matrix - matrix multiplication
    cout << "A22 = A24*A42;" << endl;
    Matrix A22 = A24*A42;
    cout << "A22 =" << endl << A22 << endl;

    // matrix - scalar operations
    cout << "A22 = A22-50.0;" << endl;
    cout << "A22 =" << endl << (A22=A22-50.0) << endl;
    cout << "A22 *= 0.5;" << endl;
    cout << "A22 =" << endl << (A22*=0.5) << endl;  // A22*=0.5 returns a const reference to A22
    cout << "A22 = A22/2.0;" << endl;
    cout << "A22 =" << endl << (A22=A22/2.0) << endl;
    cout << "A22 += 1.0;" << endl;
    // A22+=1.0 returns a const reference to A22
    cout << "A22 =" << endl << (A22+=1.0) << endl;
    cout << "A24b = A22;" << endl;
    Matrix A24b = A22;  // copy constructor is invoked
    cout << "A24b *= A24;" << endl;
    cout << "A24b =" << endl << (A24b*=A24) << endl;

    // matrix - matrix addition/subtraction
    cout << "A24b - A24 =" << endl << A24b-A24 << endl;
    cout << "A24b += A24;" << endl;
    cout << "A24b =" << endl << (A24b+=A24) << endl;

    // matrix == matrix
    cout << "A24 == A22?    " << (A24==A22) << endl;
    cout << "A24 == A24b?   " << (A24==A24b) << endl;
    cout << "A24 != A24?    " << (A24!=A24) << endl;
    cout << endl;

    // horizontal concatenation [A24,A22] or equivalently horzcat(A,B) in OCTAVE/MATLAB
    cout << "Matrix::horizontalConcatenate(A24,A22) =" << endl << Matrix::horizontalConcatenate(A24,A22) << endl;

    // vertical concatenation [A24;A22] or equivalently vertcat(A,B) in OCTAVE/MATLAB
    cout << "Matrix::verticalConcatenate(A42,A22) =" << endl << Matrix::verticalConcatenate(A42,A22) << endl;

    cout << endl << endl;


    ///////////////////////////////////////////////////////////////////////////
    //    matrix - vector operations
    ///////////////////////////////////////////////////////////////////////////
    cout << "### matrix - vector operations" << endl << endl;

    // matrix construction
    Matrix A(3,3);
    A(1,1) = 1.0;    A(1,2) = 4.0;    A(1,3) = 7.0;
    A(2,1) = 2.0;    A(2,2) = 5.0;    A(2,3) = 8.0;
    A(3,1) = 3.0;    A(3,2) = 6.0;    A(3,3) = 9.0;
    cout << "A =" << endl << A << endl;

    // matrix - vector multiplication
    cout << "v =" << endl << v << endl;
    cout << "A*v =" << endl << A*v << endl;

    // vector - matrix multiplication
    cout << "v'*A =" << endl << v*A << endl;

    cout << endl << endl;


    ///////////////////////////////////////////////////////////////////////////
    //    symmetric matrix routines
    ///////////////////////////////////////////////////////////////////////////
    cout << "### symmetric matrix routines" << endl << endl;

    // symmetric matrix assignment
    // for B = [ 1  2  4;  2  3  5;  4  5  6 ] as a symmetric matrix
    SymmetricMatrix B(3);
    B(1,1) = 1.0;
    B(2,1) = 2.0;    B(2,2) = 3.0;
    B(3,1) = 4.0;    B(3,2) = 5.0;    B(3,3) = 6.0;
    cout << "B =" << endl << B << endl;

    // convert B to a general matrix
    cout << "B.toGeneral() =" << endl << B.toGeneral() << endl;

    // submatrix of a symmetric matrix
    // for B(:,3) as a vector
    cout << "B.column(3) =" << endl << B.column(3) << endl;
    // for B(2,:) as a vector
    cout << "B.row(2) =" << endl << B.row(2) << endl;
    // for B(:,2:3) as a general matrix
    cout << "B.columns(2,3) =" << endl << B.columns(2,3) << endl;
    // for B(2:-1:1,:) as a general matrix
    cout << "B.rows(2,1) =" << endl << B.rows(2,1) << endl;
    // for B(2:3,2:-1:1) as a general matrix
    cout << "B.subMatrix(2,3,2,1) =" << endl << B.subMatrix(2,3,2,1) << endl;
    // for B(2:3,2:3) as a symmetric matrix
    cout << "B.subMatrix(2,3) =" << endl << B.subMatrix(2,3) << endl;
    // for B(2:-1:1,2:-1:1) as a symmetric matrix
    cout << "B.subMatrix(2,1) =" << endl << B.subMatrix(2,1) << endl;

    // diagonal of a symmetric matrix
    // for diag(B), which is equivalent to diag(B,0), as a vector
    cout << "B.diag() =" << endl << B.diag() << endl;
    // for diag(B,1) as a vector
    cout << "B.diag(1) =" << endl << B.diag(1) << endl;
    // for diag(B,-2) as a vector
    cout << "B.diag(-2) =" << endl << B.diag(-2) << endl;

    // transpose of a symmetric matrix
    // for B' as a symmetric matrix
    cout << "B.transpose() =" << endl << B.transpose() << endl;

    // norms of a symmetric matrix
    cout << "||B||_F^2 = " << Basics::square(B.normFrobenius()) << endl;
    cout << "||B||_1   = " << B.norm1() << endl;
    cout << "||B||_inf = " << B.normInfinity() << endl << endl;

    // symmetric permutation of rows and columns of a symmetric matrix
    cout << "perm3 = [ 2, 3, 1 ]" << endl << endl;
    // for B(perm3,perm3) as a symmetric matrix
    SymmetricMatrix C = B.permuteRowsAndColumns(perm3);
    cout << "C = B.permuteRowsAndColumns(perm3);" << endl;
    cout << "C = " << endl << C << endl;
    // cout << "perm3x = [ 1, 2, 0 ]" << endl;
    // cout << "B.permuteRowsAndColumns(perm3x,false) =" << endl << B.permuteRowsAndColumns(perm3x,false) << endl;  // same result as above

    // symmetric matrix == symmetric matrix
    cout << "B == C?    " << (B==C) << endl;
    cout << "B != C?    " << (B!=C) << endl;
    cout << endl;

    // symmetric matrix - symmetric matrix addition/subtraction
    cout << "B+C =" << endl << B+C << endl;
    cout << "B-C =" << endl << B-C << endl;

    // symmetric matrix - symmetric matrix multiplication
    cout << "B*C =" << endl << B*C << endl;

    cout << endl << endl;


    ///////////////////////////////////////////////////////////////////////////
    //    symmetric matrix - vector operations
    ///////////////////////////////////////////////////////////////////////////
    cout << "### symmetric matrix - vector operations" << endl << endl;

    // symmetric matrix - vector multiplication
    cout << "v =" << endl << v << endl;
    cout << "B*v =" << endl << B*v << endl;

    // vector - symmetric matrix multiplication
    cout << "v'*B =" << endl << v*B << endl;

    cout << endl << endl;


    ///////////////////////////////////////////////////////////////////////////
    //    matrix - symmetric matrix operations
    ///////////////////////////////////////////////////////////////////////////
    cout << "### matrix - symmetric matrix operations" << endl << endl;

    // matrix - symmetric matrix addition/subtraction
    cout << "A =" << endl << A << endl;
    cout << "A+B =" << endl << A+B << endl;
    cout << "A -= B;" << endl;
    cout << "A =" << endl << (A-=B) << endl;

    // symmetric matrix - matrix addition/subtraction
    cout << "B+A =" << endl << B+A << endl;
    cout << "B-A =" << endl << B-A << endl;

    // matrix assignment and transpose
    // for A23 = [ 1  4  7;  2  4  8 ]
    Matrix A23(2,3);
    A23(1,1) = 1.0;    A23(1,2) = 4.0;    A23(1,3) = 7.0;
    A23(2,1) = 2.0;    A23(2,2) = 6.0;    A23(2,3) = 8.0;
    cout << "A23 =" << endl << A23 << endl;
    cout << "A32 = A23.transpose();" << endl;
    Matrix A32 = A23.transpose();
    cout << "A32 =" << endl << A32 << endl;

    // symmetric matrix - matrix multiplication
    cout << "B*A32 =" << endl << B*A32 << endl;

    // matrix - symmetric matrix multiplication
    cout << "A23*B =" << endl << (A23*B) << endl;
    cout << "A23 *= B;" << endl;
    cout << "A23 =" << endl << (A23*=B) << endl;

    cout << endl << endl;


    ///////////////////////////////////////////////////////////////////////////
    //    sparse matrix routines
    ///////////////////////////////////////////////////////////////////////////
    cout << "### sparse matrix routines" << endl << endl;

    // speye function and scaler multiplication
    SparseMatrix S22 = SparseMatrix::eye(2);
    S22 *= 3.0;
    cout << "S22 =" << endl << S22 << endl;

    // read elements via a const reference
    const SparseMatrix &S22c = S22;
    cout << "S22(1,2) == " << S22c(1,2) << endl;
    cout << "S22(2,2) == " << S22c(2,2) << endl << endl;

    // write elements
    SparseMatrix S32 = SparseMatrix::eye(3,2);
    S32 /= 2.0;
    cout << "S32 =" << endl << S32 << endl;
    cout << "S32(2,2) *= 2.0;" << endl;
    S32(2,2) *= 2.0;
    cout << "S32(3,1) = 5.0;" << endl;
    S32(3,1) = 5.0;
    cout << endl << "S32 =" << endl << S32 << endl;

    // full(S32)
    cout << "S32.toFull()" << endl << S32.toFull() << endl;

    // Frobenius norm of a sparse matrix
    cout << "||S32||_F^2 = " << Basics::square(S32.normFrobenius()) << endl;
    cout << endl;

    // diagonal of a sparse matrix (stored as a vector)
    // full(diag(S32))
    cout << "S32.diag() = " << S32.diag() << endl;
    // full(diag(S32,1))
    cout << "S32.diag(1) = " << S32.diag(1) << endl;
    // full(diag(S32,-2))
    cout << "S32.diag(-2) = " << S32.diag(-2) << endl;

    // write/read elements
    cout << "S32(3,2) = -1.0;" << endl;
    S32(3,2) = -1.0;
    const SparseMatrix& S32c = S32;
    cout << "display: S32(3,2) = " <<  S32c(3,2) << ", S32(2,2) = " << S32c(2,2) << endl;
    cout << endl;

    // sparse matrix transpose
    cout << "S32 =" << endl << S32 << endl;
    cout << "S32.transpose() =" << endl << S32.transpose() << endl;

    // permutation of rows / columns of a sparse matrix
    cout << "perm2 = [ 2, 1 ]" << endl;
    // cout << "perm2x = [ 1, 0 ]" << endl;
    cout << "perm3 = [ 2, 3, 1 ]" << endl;
    // cout << "perm3x = [ 1, 2, 0 ]" << endl;
    cout << endl;
    cout << "S32.permuteColumns(perm2) =" << endl << S32.permuteColumns(perm2) << endl;
    // cout << "S32.permuteColumns(perm2x, false) =" << endl << S32.permuteColumns(perm2x, false) << endl;  // same result as above
    cout << "S32.permuteRows(perm3) =" << endl << S32.permuteRows(perm3) << endl;
    // cout << "S32.permuteRows(perm3x, false) =" << endl << S32.permuteRows(perm3x, false) << endl;  // same result as above
    cout << "S32.permuteRowsAndColumns(perm3, perm2) =" << endl << S32.permuteRowsAndColumns(perm3, perm2) << endl;
    // cout << "S32.permuteRowsAndColumns(perm3x, perm2x, false) =" << endl << S32.permuteRowsAndColumns(perm3x, perm2x, false) << endl;  // same result as above

    // sparse matrix assignment
    SparseMatrix T32(3,2);
    T32(2,1) = 0.2;
    T32(1,2) = 0.5;
    T32(2,2) = 2.1;
    cout << "T32 =" << endl;
    cout << T32 << endl;

    // sparse matrix - sparse matrix addition / subtraction
    cout << "U32 = S32 + T32;" << endl;
    SparseMatrix U32 = S32 + T32;
    cout << "U32 =" << endl << U32 << endl;
    cout << "T32 -= S32;" << endl;
    cout << "T32 =" << endl << (T32 -= S32) << endl;

    // submatrix of a sparse matrix
    // for U32(:,2:-1:1) as a sparse matrix
    cout << "U32.columns(2,1) =" << endl << U32.columns(2,1) << endl;
    // for U32(3:-1:2,:) as a sparse matrix
    cout << "U32.rows(3,2) =" << endl << U32.rows(3,2) << endl;
    // for U32(3:-1:2,2:-1:1) as a sparse matrix
    cout << "U32.subMatrix(3,2,2,1) =" << endl << U32.subMatrix(3,2,2,1) << endl;

    // sparse matrix == sparse matrix
    cout << "T32 == S32?    " << (T32==S32) << endl;
    cout << "T32 != S32?    " << (T32!=S32) << endl;
    cout << "T32 == T32?    " << (T32==T32) << endl;
    cout << "T32 != T32?    " << (T32!=T32) << endl << endl;

    // form a 3-by-3 tridiagonal matrix as a sparse matrix
    // T = sparse(diag([2 3],-1) + diag([3 4 5]) + diag([1 2],1))
    Real st[] = { 1.0, 2.0, 3.0, 4.0, 5.0 };
    SparseMatrix T = SparseMatrix::tridiagonalMatrix(3, st+1, st+2, st);
    cout << "T =" << endl << T << endl;

    cout << endl << endl;


    ///////////////////////////////////////////////////////////////////////////
    //    sparse matrix - vector operations
    ///////////////////////////////////////////////////////////////////////////
    cout << "### sparse matrix - vector operations" << endl << endl;

    // sparse matrix - vector (and vector - sparse matrix) multiplication
    v.resize(2);
    v(1) = 1.0;  v(2) = 2.0;
    cout << "v =" << endl << v << endl;
    w = S32*v;
    cout << "w = S32*v;" << endl;
    cout << "w =" << endl << w << endl;
    cout << "w *= S32;" << endl;
    w *= S32;
    cout << "w =" << endl << w << endl;

    cout << endl << endl;


    ///////////////////////////////////////////////////////////////////////////
    //    matrix - sparse matrix operations
    ///////////////////////////////////////////////////////////////////////////
    cout << "### matrix - sparse matrix operations" << endl << endl;

    // sparse matrix - matrix multiplication
    cout << "S32 =" << endl << S32 << endl;
    cout << "A23 =" << endl << A23 << endl;
    cout << "S32*A23 =" << endl << S32*A23 << endl;

    // matrix - sparse matrix multiplication
    cout << "A23*S32 =" << endl << A23*S32 << endl;

    return 0;
}
