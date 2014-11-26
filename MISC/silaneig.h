#ifndef SILANEIG_H
#define SILANEIG_H

#include "laneig.h"
#include "umfpack.h"

#ifdef USE_NAMESPACE
using namespace MATKIT;
namespace MATKIT {
#endif

class Vector;
class Matrix;
class SymmetricMatrix;
class SparseMatrix;

#ifdef USE_NAMESPACE
}  // end of namespace MATKIT
#endif


class ShiftAndInvertInterface {
private:
    // works for sparse matrices only
    static const SparseMatrix *Sptr;   
    static SparseMatrix S2;   

    // variables for UMFPACK parameter setting
    static double UMFinfo[UMFPACK_INFO];
    static double UMFcontrol[UMFPACK_CONTROL];

    // variables for the matrix used by UMFPACK
    static UF_long *Ap, *Ai;

    // variables for sparse LU of the shifted matrix by UMFPACK
    static void *UMFsymbolic;
    static void *UMFnumeric; 

public:
    // routines for shift-and-invert Lanczos
    static void setLUbyUMF(const SparseMatrix &S, Real shift);
    static Vector ShiftAndInvertSparseMatrixVectorProduct(const Vector &v);
    static void freeUMF();
};

struct ShiftAndInvertLanczosOptions: LanczosOptions {
    Real eigSortBeforeShiftAndInvert;
    // -2 for "don't sort"
    // -1 for default setting (0 for increasing; the corresponding LanczosOptions::eigSort is 2 for "LA"; 1 for "SA"; 3 for "BE")
    //        note that LanczosOptions::eigSort indicates the sort order after shift and invert
    //  otherwise eigSortBeforeShiftAndInvert should be >= 0, in which case
    //        the first digit is 0 (or 1) for the eigenvalues in the lower end in increasing (or decreasing) order
    //        the second digit is 0 (or 1) for the eigenvalues in the higher end in increasing (or decreasing) order
    //        the third digit is 0 (or 1) for the eigenvalues in the lower end (or higher end) recorded first
    //  note: when eigPart[] is "LA", the first and the third digits do not matter
    //        when eigPart[] is "SA", the second and the third digits do not matter
    //        when eigPart[] is "BE", set eigSort=0 (or 7) for eigenvalues sorted in increasing (or decreasing) order

    Real eigLowerBound, eigUpperBound;

    // default constructor
    ShiftAndInvertLanczosOptions() {
        eigSortBeforeShiftAndInvert = -1;    // the default sort order
        eigLowerBound = -Basics::infinity;   // means no bound
        eigUpperBound = Basics::infinity;    // means no bound
        defaultMinIterFactor = 2;
        defaultMaxIterFactor = 30;
    }
};

struct ShiftAndInvertLanczosInfo: LanczosInfo {
    Real LuCpuTime;
};

// Shift-and-invert Lanczos eigensolver
ShiftAndInvertLanczosInfo ShiftAndInvertLanczosEigenSolver(Vector &eigVal, Matrix &eigVec,
                          const SparseMatrix &S, Real shift, mkIndex neigWanted, const char eigPart[], ShiftAndInvertLanczosOptions &opts);
ShiftAndInvertLanczosInfo ShiftAndInvertLanczosEigenSolver(Vector &eigVal, Matrix &eigVec,
                          const SparseMatrix &S, Real shift, mkIndex neigWanted, const char eigPart[]);  // default opts will be used

#endif
