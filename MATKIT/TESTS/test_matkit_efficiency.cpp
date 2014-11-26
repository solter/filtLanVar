// this driver is for testing the efficiency of some MATKIT functions

#include <iostream>  // for cout, cerr, endl, etc. (under namespace std)
#include <time.h>    // for time_t, time, clock_t, clock, CLOCKS_PER_SEC, etc.
#include "matkit.h"

using std::cout;
using std::endl;

#ifdef USE_NAMESPACE
using namespace MATKIT;
#endif


int main(int argc, char *argv[]) {
    ///////////////////////////////////////////////////////////////////////////
    //    efficiency test
    ///////////////////////////////////////////////////////////////////////////
    cout << "### efficiency test" << endl << endl;

    // construct a sparse matrix for tests
    mkIndex nr = 20000, nc = 30000;
    int bnd = 50;
    SparseMatrix S(nr, nc);
    for (int j=1; j<=(int)nc; j++) {
        for (int i=Basics::max(j-bnd,1); i<=Basics::min(j+bnd,(int)nr); i++) {
            S(i,j) = 1.0/(1.0+(Real)Basics::abs(i-j));
        }
    }
    cout << "S is of size " << S.Nrows() << "-by-" << S.Ncols() << " with nnz=" << S.Nnz() << endl;
    cout << "||S||_1   = " << S.norm1() << endl;
    cout << "||S||_F   = " << S.normFrobenius() << endl;
    cout << "||S||_inf = " << S.normInfinity() << endl << endl;

    // sparse matrix - vector multiplication
    cout << "v = Vector::ones(" << S.Ncols() << ")" << endl;
    Vector v = Vector::ones(S.Ncols()), w;
    cout << "w=S*v time used (1000 times): ";
    clock_t start = clock();
    for (int i=0; i<1000; i++)
        w = S*v;
    cout << (Real)(clock()-start)/(Real)CLOCKS_PER_SEC << " seconds" << endl;
    cout << "||w||_1   = " << w.norm1() << endl;
    cout << "||w||_2   = " << w.norm2() << endl;
    cout << "||w||_inf = " << w.normInfinity() << endl << endl;

    cout << "v = Vector::ones(" << S.Nrows() << ")" << endl;
    v = Vector::ones(S.Nrows());
    Vector w2;
    cout << "w2=v'*S time used (1000 times): ";
    start = clock();
    for (int i=0; i<1000; i++)
        w2 = v*S;
    cout << (Real)(clock()-start)/(Real)CLOCKS_PER_SEC << " seconds" << endl;
    cout << "||w2||_1   = " << w2.norm1() << endl;
    cout << "||w2||_2   = " << w2.norm2() << endl;
    cout << "||w2||_inf = " << w2.normInfinity() << endl << endl;

    // dense matrix transpose
    mkIndex dim = 10000;
    Matrix M = Matrix::random(dim, dim);
    cout << "time used to transpose a " << dim << "-by-" << dim << " dense matrix: ";
    start = clock();
    M.transpose();
    cout << (Real)(clock()-start)/(Real)CLOCKS_PER_SEC << " seconds" << endl << endl;

    // dense matrix - vector multiplication
    v = Vector::random(dim);
    start = clock();
    cout << "dense matrix - vector multiplication (n=" << dim << ") time used: ";
    v = M*v;
    cout << (Real)(clock()-start)/(Real)CLOCKS_PER_SEC << " seconds" << endl;

    cout << "vector - dense matrix multiplication (n=" << dim << ") time used: ";
    start = clock();
    v = v*M;
    cout << (Real)(clock()-start)/(Real)CLOCKS_PER_SEC << " seconds" << endl;

    SymmetricMatrix M2 = SymmetricMatrix::random(dim);
    cout << "symmetric matrix - vector multiplication (n=" << dim << ") time used: ";
    start = clock();
    v = M2*v;
    cout << (Real)(clock()-start)/(Real)CLOCKS_PER_SEC << " seconds" << endl << endl;

    // matrix - matrix multiplication, matrix - symmetric matrix multiplication
    dim = 2000;
    M = Matrix::random(dim, dim);
    cout << "dense matrix - dense matrix multiplication (n=" << dim << ") time used: ";
    start = clock();
    Matrix N = M*M;
    cout << (Real)(clock()-start)/(Real)CLOCKS_PER_SEC << " seconds" << endl;
    M2 = SymmetricMatrix::random(dim);
    cout << "dense matrix - symmetric matrix multiplication (n=" << dim << ") time used: ";
    start = clock();
    N = M*M2;
    cout << (Real)(clock()-start)/(Real)CLOCKS_PER_SEC << " seconds" << endl;
    cout << "symmetric matrix - dense matrix multiplication (n=" << dim << ") time used: ";
    start = clock();
    N = M2*M;
    cout << (Real)(clock()-start)/(Real)CLOCKS_PER_SEC << " seconds" << endl;

    return 0;
}
