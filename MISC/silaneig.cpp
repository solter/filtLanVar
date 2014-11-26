#include <stdlib.h>  // for exit, atof, atoi, atol, srand, rand, etc.
#include <string.h>  // for strcmp
#include <time.h>    // for time_t, time, clock_t, clock, CLOCKS_PER_SEC, etc.
#include <iostream>  // for cout, cerr, endl, etc. (under namespace std)

#include "matkit.h"
#include "laneig.h"
#include "silaneig.h"
#include "umfpack.h"

using std::cerr;
using std::endl;


////////////////////////////////////////////////////////////////////////////////

// static variables for ShiftAndInvertInterface
const SparseMatrix *ShiftAndInvertInterface::Sptr = NULL;
SparseMatrix ShiftAndInvertInterface::S2;

// variables for UMFPACK parameter setting
double ShiftAndInvertInterface::UMFinfo[UMFPACK_INFO];
double ShiftAndInvertInterface::UMFcontrol[UMFPACK_CONTROL];

// variables for the matrix used by UMFPACK
UF_long *ShiftAndInvertInterface::Ap = NULL;
UF_long *ShiftAndInvertInterface::Ai = NULL;

// variables for the sparse LU of the shifted matrix
void *ShiftAndInvertInterface::UMFsymbolic = NULL;
void *ShiftAndInvertInterface::UMFnumeric = NULL;

// static member functions of ShiftAndInvertInterface
void ShiftAndInvertInterface::setLUbyUMF(const SparseMatrix &S, Real shift) {
    // get the shifted matrix
    mkIndex n = S.Nrows();
    if (n != S.Ncols()) {
        cerr << "ShiftAndInvertInterface::setLUbyUMF(const SparseMatrix &, Real): the given matrix is not square!" << endl;
        exit(1);
    }
    S2 = S - shift*SparseMatrix::eye(n, n);

    // get the default UMF control parameters
    umfpack_dl_defaults(UMFcontrol);

    // convert the indices, which describe the sparsity pattern, to be UF_long for UMFPACK to read
    mkIndex nnz = S2.Nnz();
    Ap = new UF_long[n+1];
    Ai = new UF_long[nnz];
    UF_long *ap = Ap;
    UF_long *ai = Ai;
    const mkIndex *sp = S2.ColPointer();
    const mkIndex *si = S2.RowIndex();
    mkIndex j = n+1;
    while (j--)
        *ap++ = *sp++;
    j = nnz;
    while (j--)
        *ai++ = *si++;

    // compute LU of the shifted spmatrix S2 by UMFPACK
    UF_long status = umfpack_dl_symbolic((UF_long)n, (UF_long)n, Ap, Ai, S2.Store(), &UMFsymbolic, UMFcontrol, UMFinfo);
    if (status < 0) {
        umfpack_dl_report_info(UMFcontrol, UMFinfo);
        umfpack_dl_report_status(UMFcontrol, status);
        cerr << "ShiftAndInvertInterface::setLUbyUMF(const SparseMatrix &, Real): umfpack_dl_symbolic failed!" << endl;
        exit(1);
    }
    status = umfpack_dl_numeric(Ap, Ai, S2.Store(), UMFsymbolic, &UMFnumeric, UMFcontrol, UMFinfo);
    if (status < 0) {
        umfpack_dl_report_info(UMFcontrol, UMFinfo);
        umfpack_dl_report_status(UMFcontrol, status);
        cerr << "ShiftAndInvertInterface::setLUbyUMF(const SparseMatrix &, Real): umfpack_dl_numeric failed!" << endl;
        exit(1);
    }
}

Vector ShiftAndInvertInterface::ShiftAndInvertSparseMatrixVectorProduct(const Vector &v) {
    mkIndex n = v.Length();
    if (n != S2.Ncols()) {
        cerr << "ShiftAndInvertInterface::ShiftAndInvertSparseMatrixVectorProduct(const Vector &): inner matrix dimensions must agree!" << endl;
        exit(1);
    }
    Real *b = new Real[n];
    UF_long status = umfpack_dl_solve(UMFPACK_A, Ap, Ai, S2.Store(), b, v.Store(), UMFnumeric, UMFcontrol, UMFinfo);
    if (status < 0) {
        cerr << "ShiftAndInvertInterface::ShiftAndInvertSparseMatrixVectorProduct(const Vector &): umfpack_dl_solve failed!" << endl;
        exit(1);
    }
    return Vector(n, b);
}

void ShiftAndInvertInterface::freeUMF() {
    umfpack_dl_free_symbolic(&UMFsymbolic);
    umfpack_dl_free_numeric(&UMFnumeric);
    delete [] Ai;
    Ai = NULL;
    delete [] Ap;
    Ap = NULL;
}

////////////////////////////////////////////////////////////////////////////////

ShiftAndInvertLanczosInfo ShiftAndInvertLanczosEigenSolver(Vector &eigVal, Matrix &eigVec, const SparseMatrix &S, Real shift, mkIndex neigWanted, const char eigPart[]) {
    ShiftAndInvertLanczosOptions opts;
    return ShiftAndInvertLanczosEigenSolver(eigVal, eigVec, S, shift, neigWanted, eigPart, opts);
}
ShiftAndInvertLanczosInfo ShiftAndInvertLanczosEigenSolver(Vector &eigVal, Matrix &eigVec, const SparseMatrix &S, Real shift, mkIndex neigWanted, const char eigPart[], ShiftAndInvertLanczosOptions &opts) {
    ShiftAndInvertLanczosInfo outputInfo;

    // deal with exceptions
    if (!strcmp(eigPart,"LA") || !strcmp(eigPart,"la")) {
        if (shift > opts.eigUpperBound) {
            cerr << "ShiftAndInvertLanczosEigenSolver(Vector &, Matrix &, const SparseMatrix &, Real, mkIndex, char, ShiftAndInvertLanczosOptions &): the shift value " << shift << " is greater than the upper bound " << opts.eigUpperBound << " of desired eigenvalues (eigen-part \"LA\")!" << endl;
            exit(1);
        }
    }
    else if (!strcmp(eigPart,"SA") || !strcmp(eigPart,"sa")) {
        if (shift < opts.eigLowerBound) {
            cerr << "ShiftAndInvertLanczosEigenSolver(Vector &, Matrix &, const SparseMatrix &, Real, mkIndex, char, ShiftAndInvertLanczosOptions &): the shift value " << shift << " is less than the lower bound " << opts.eigLowerBound << " of desired eigenvalues (eigen-part \"SA\")!" << endl;
            exit(1);
        }
    }
    else if (!strcmp(eigPart,"BE") || !strcmp(eigPart,"be")) {
        if (opts.eigLowerBound >= opts.eigUpperBound) {
            cerr << "ShiftAndInvertLanczosEigenSolver(Vector &, Matrix &, const SparseMatrix &, Real, mkIndex, char, ShiftAndInvertLanczosOptions &): the lower bound (" << opts.eigLowerBound << ") should be less than the upper bound (" << opts.eigUpperBound << ") (eigen-part \"BE\")!" << endl;
            exit(1);
        }
        if (shift > opts.eigUpperBound) {
            cerr << "ShiftAndInvertLanczosEigenSolver(Vector &, Matrix &, const SparseMatrix &, Real, mkIndex, char, ShiftAndInvertLanczosOptions &): the shift value " << shift << " is greater than the upper bound " << opts.eigUpperBound << " of desired eigenvalues (eigen-part \"BE\")!" << endl;
            exit(1);
        }
        if (shift < opts.eigLowerBound) {
            cerr << "ShiftAndInvertLanczosEigenSolver(Vector &, Matrix &, const SparseMatrix &, Real, mkIndex, char, ShiftAndInvertLanczosOptions &): the shift value " << shift << " is less than the lower bound " << opts.eigLowerBound << " of desired eigenvalues (eigen-part \"BE\")!" << endl;
            exit(1);
        }
    }
    else {
        cerr << "ShiftAndInvertLanczosEigenSolver(Vector &, Matrix &, const SparseMatrix &, Real, mkIndex, char, ShiftAndInvertLanczosOptions &): the eigenvalue part of spectrum \"" << eigPart << "\" is not supported!" << endl;
        exit(1);
    }

    // set sort order of eigenvalues
    if (opts.eigSortBeforeShiftAndInvert <= -2) {
        opts.eigSort = -2;  // means "don't sort"
    }
    else if (opts.eigSortBeforeShiftAndInvert == -1) {
        // the default case
        opts.eigSortBeforeShiftAndInvert = 0;  // sort in ascending order
        if (!strcmp(eigPart,"LA") || !strcmp(eigPart,"la"))
            opts.eigSort = 2;  // the transformed eigenvalues are decreasing; the wanted eigenvalues are increasing
        else if (!strcmp(eigPart,"SA") || !strcmp(eigPart,"sa"))
            opts.eigSort = 1;  // the transformed eigenvalues are decreasing; the wanted eigenvalues are increasing
        else if (!strcmp(eigPart,"BE") || !strcmp(eigPart,"be"))
            opts.eigSort = 3;  // the transformed eigenvalues in the lower (or higher) end are decreasing;
                               // the transformed eigenvalues are decreasing; the wanted eigenvalues are increasing
    }
    else {
        // the general case
        opts.eigSort = opts.eigSortBeforeShiftAndInvert;
        // reverse the order of transformed eigenvalues in the lower end
        if (opts.eigSort%2 == 1)
            opts.eigSort -= 1;
        else
            opts.eigSort += 1;
        // reverse the order of transformed eigenvalues in the higher end
        if ((opts.eigSort/2)%2 == 1)
            opts.eigSort -= 2;
        else
            opts.eigSort += 2;
        // in either end, the wanted eigenvalues and the transformed eigenvalues are in reverse order
        // so we will get the wanted eigenvalues in the right order
    }

    // set the eigen cuts (after shift-and-invert) according to the eigen bounds (before shift-and-invert)
    if (opts.eigLowerBound > -Basics::infinity) {
        // effective if !strcmp(eigPart,"SA") || !strcmp(eigPart,"sa") || !strcmp(eigPart,"BE") || !strcmp(eigPart,"BE")
        opts.eigLowCut = 1.0 / (opts.eigLowerBound-shift);   // transformed eigen cut
    }
    if (opts.eigUpperBound < Basics::infinity) {
        // effective if !strcmp(eigPart,"LA") || !strcmp(eigPart,"la") || !strcmp(eigPart,"BE") || !strcmp(eigPart,"BE")
        opts.eigHighCut = 1.0 / (opts.eigUpperBound-shift);  // transformed eigen cut
    }

    // compute the LU factorization of the shifted matrix by UMFPACK
    clock_t start = clock();
    ShiftAndInvertInterface::setLUbyUMF(S, shift);
    outputInfo.LuCpuTime = (Real)(clock()-start)/(Real)CLOCKS_PER_SEC;

    // compute the eigenvectors
    LanczosInfo &shadowInfo = outputInfo;
    shadowInfo = LanczosEigenSolver(eigVal, eigVec, ShiftAndInvertInterface::ShiftAndInvertSparseMatrixVectorProduct, S.Nrows(), neigWanted, eigPart, opts);

    // compute the eigenvalues
    mkIndex neig = eigVal.Length();
    for (mkIndex i=1; i<=neig; i++)
        eigVal(i) = shift + 1.0/eigVal(i);  // transform the eigenvalues back

    // free memory used by UMFPACK
    ShiftAndInvertInterface::freeUMF();

    return outputInfo;
}
