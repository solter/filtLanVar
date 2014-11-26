//============================================================================
// The filtlan mex driver for OCTAVE/MATLAB users
//
// Reference:
// "A Filtered Lanczos Procedure for Extreme and Interior Eigenvalue Problems"
// H.-r. Fang and Y. Saad, University of Minnesota Technical Report, 2011.
//============================================================================

#include <stdio.h>
#include <string.h>
#include "matkit.h"
#include "filtlan.h"
#include "mex.h"

void printUsage();
void parseFilteredLanczosOptions(FilteredLanczosOptions &opts, const mxArray *mxOpts);
mxArray *convertFilteredLanczosInfo(const FilteredLanczosInfo &info);


// OCTAVE/MATLAB: [V, D, info] = filtlan(A, [a,b], opts);
//                [V, D]       = filtlan(A, [a,b], opts);
//                    D        = filtlan(A, [a,b], opts);
// OCTAVE/MATLAB: [V, D, info] = filtlan(A, [a,b], polydeg, opts);
//                [V, D]       = filtlan(A, [a,b], polydeg, opts);
//                    D        = filtlan(A, [a,b], polydeg, opts);
// OCTAVE/MATLAB: [V, D, info] = filtlan(A, [a,b], polydeg, basedeg, opts);
//                [V, D]       = filtlan(A, [a,b], polydeg, basedeg, opts);
//                    D        = filtlan(A, [a,b], polydeg, basedeg, opts);
// ======
// input:
// A is a sparse matrix whose eigenvalues are of interest.
// [a,b] is the interval in which the eigenvalues are requested
//    if a==-inf, then the requested eigenvalues are in the lower end of the spectrum
//    if b== inf, then the requested eigenvalues are in the upper end of the spectrum
//    if [a,b] is properly inside the spectrum, the requested eigenvalues are interior
// polydeg is the polynomial degree (default 10, which should be increased if the requested
//    eigenvalues are interior especially in a narrow interval relative to the range of spectrum)
// basedeg is the base filter degree (default 10); degrees from 5 to 15 are usually good
// opts is a collection of filtered Lanczos options, which is optional, and
// by default it is from the default constructor of FilteredLanczosOptions
// ======
// output:
// V is a dense matrix whose columns are the computed eigenvectors of A
// D is a vector or a sparse diagonal matrix formed by the computed eigenvalues of A
// info gives some information about the computation
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    if (nrhs == 0) {
        printUsage();
        return;
    }
    if (nrhs < 2)
        mexErrMsgTxt("filtlan: at least 2 input arguments should be provided!");
    else if (nrhs > 5)
        mexErrMsgTxt("filtlan: at most 5 input arguments should be provided!");
    if (nlhs > 3)
        mexErrMsgTxt("filtlan: at most 3 output arguments are generated!");

    // OCTAVE/MATLAB: [nrows, ncols] = size(A);
    mkIndex nrows = (mkIndex)mxGetM(prhs[0]);
    mkIndex ncols = (mkIndex)mxGetN(prhs[0]);
    if (nrows != ncols)
        mexErrMsgTxt("filtlan: input matrix must be square!");

    if (!mxIsSparse(prhs[0]))
        mexErrMsgTxt("filtlan: input as a dense matrix is not (yet) supported!");

    SparseMatrix A = SparseMatrix::mxArray2SparseMatrix(prhs[0]);
    // the routine mxArray2SparseMatrix is declared in spmatrix.h of MATKIT
    // assume that the input matrix A is symmetric

    // get the interval of eigenvalues requested
    Vector frame(2);
    if (!mxIsNumeric(prhs[1]) || mxGetM(prhs[1])*mxGetN(prhs[1])!=2)
        mexErrMsgTxt("filtlan: the 2nd input argument must contain exactly two numbers specifying the interval of the eigenvalues!");
    frame(1) = mxArray_to_double(prhs[1]);
    frame(2) = mxArray_to_double(prhs[1], 1);

    // get the requested number of eigenvalues, if it is specified as the 3rd input argument
    // also set options, if it is specified, from the 3rd or the 4th input argument
    FilteredLanczosOptions opts;  // the default is set by the constructor
    mkIndex polyDeg = 10, baseDeg = 10;
    if (nrhs >= 3) {
        if (mxIsNumeric(prhs[2])) {
            if (mxGetM(prhs[2])!=1 || mxGetN(prhs[2])!=1)
                mexErrMsgTxt("filtlan: the 3rd input argument should be either a single number as the polynomial degree, or a structure array for filtered Lanczos options!");
            double deg = mxArray_to_double(prhs[2]);
            if (deg < 0.0)
                mexErrMsgTxt("filtlan: the 3rd input argument, when it is the polynomial degree, should be a nonnegative integer!");
            polyDeg = (mkIndex)deg;
        }
        else if (!mxIsStruct(prhs[2])) {
            mexErrMsgTxt("filtlan: the 3rd input argument should be either a single number as the polynomial degree, or a structure array for filtered Lanczos options!");
        }
    }
    if (nrhs >= 4) {
        if (mxIsNumeric(prhs[3])) {
            if (mxGetM(prhs[3])!=1 || mxGetN(prhs[3])!=1)
                mexErrMsgTxt("filtlan: the 4th input argument should be either a single number as the base filter degree, or a structure array for filtered Lanczos options!");
            double deg = mxArray_to_double(prhs[3]);
            if (deg < 1.0)
                mexErrMsgTxt("filtlan: the 4th input argument, when it is the base filter degree, should be a positive integer!");
            baseDeg = (mkIndex)deg;
        }
        else if (!mxIsStruct(prhs[3])) {
            mexErrMsgTxt("filtlan: the 4th input argument should be either a single number as the base filter degree, or a structure array for filtered Lanczos options!");
        }
    }
    if (mxIsStruct(prhs[nrhs-1])) {
        parseFilteredLanczosOptions(opts, prhs[nrhs-1]);
    }

    // find the eigenvalues/eigenvectors requested
    Vector lambda;
    Matrix V;
    FilteredLanczosInfo outputInfo = FilteredLanczosEigenSolver(lambda, V, A, frame, polyDeg, baseDeg, opts);

    // set the output
    if (nlhs <= 1) {
        plhs[0] = Vector::Vector2mxArray(lambda);
    }
    else {  // nlhs >= 2
        plhs[0] = Matrix::Matrix2mxArray(V);  // 1st output argument
        plhs[1] = SparseMatrix::SparseMatrix2mxArray(lambda.spdiag());  // 2nd output argument
        if (nlhs == 3) {
            // now set the 3rd output argument
            plhs[2] = convertFilteredLanczosInfo(outputInfo);
        }
    }
}

void parseFilteredLanczosOptions(FilteredLanczosOptions &opts, const mxArray *mxOpts) {
    // get input arguments
    int numFields = mxGetNumberOfFields(mxOpts);
    if (mxGetNumberOfElements(mxOpts) != 1) {
        mexErrMsgTxt("filtlan options parser: each field should have one and only one element!");
    }

    // parse the options
    for (int ifield=0; ifield<numFields; ifield++) {
        size_t nrows, ncols;
        mxArray *fieldData = mxGetFieldByNumber(mxOpts, 0, ifield);
        if ((nrows=mxGetM(fieldData)) == 0 || (ncols=mxGetN(fieldData)) == 0)
            continue;  // empty field, so ignore
        const char *fieldName = mxGetFieldNameByNumber(mxOpts, ifield);
        if (!strcmp(fieldName, "v0")) {
            opts.v0 = Vector::mxArray2Vector(fieldData);
            continue;
        }
        if (!strcmp(fieldName, "intervalWeights")) {
            opts.intervalOpts.intervalWeights = Vector::mxArray2Vector(fieldData);
            continue;
        }
        double val = mxArray_to_double(fieldData);
        if (!strcmp(fieldName, "rseed")) {
            if (nrows!=1 || ncols!=1 || val<0.0)
                mexErrMsgTxt("filtlan options error: the random seed, rseed, should be a single nonnegative integer!");
            srand((unsigned)val);
        }
        else if (!strcmp(fieldName, "eigrangeiter")) {
            if (nrows!=1 || ncols!=1 || val<0.0)
                mexErrMsgTxt("filtlan options error: the number of Lanczos iterations to determine the eigen-range, eigrangeiter, should be a single nonnegative number!");
            opts.numIterForEigenRange = (mkIndex)val;
        }
        else if (!strcmp(fieldName, "maxiter")) {
            if (nrows!=1 || ncols!=1 || val<0.0)
                mexErrMsgTxt("filtlan options error: maximum number of iterations, maxiter, should be a single nonnegative number!");
            opts.maxIter = (mkIndex)val;
        }
        else if (!strcmp(fieldName, "miniter")) {
            if (nrows!=1 || ncols!=1 || val<0.0)
                mexErrMsgTxt("filtlan options error: minimum number of iterations, miniter, should be a single nonnegative number!");
            opts.minIter = (mkIndex)val;
        }
        else if (!strcmp(fieldName, "extraiter")) {
            if (nrows!=1 || ncols!=1 || val<0.0)
                mexErrMsgTxt("filtlan options error: extra number of iterations (after deemed convergence), extraiter, should be a single nonnegative number!");
            opts.extraIter = (mkIndex)val;
        }
        else if (!strcmp(fieldName, "stride")) {
            if (nrows!=1 || ncols!=1 || val<0.0)
                mexErrMsgTxt("filtlan options error: the convergence to be checked every 'stride' iterations, the stride should be a single nonnegative number!");
            opts.stride = (mkIndex)val;
        }
        else if (!strcmp(fieldName, "tol")) {
            if (nrows!=1 || ncols!=1)
                mexErrMsgTxt("filtlan options error: the tolerance for checking convergence of eigenvalues, tol, should be a single number!");
            opts.tol = (Real)val;
        }
        else if (!strcmp(fieldName, "disp")) {
            if (nrows!=1 || ncols!=1 || val<0.0)
                mexErrMsgTxt("filtlan options error: the diagnostic information display level, disp, should be a single nonnegative number!");
            opts.disp = (int)val;
        }
        else if (!strcmp(fieldName, "expand")) {
            if (nrows!=1 || ncols!=1)
                mexErrMsgTxt("filtlan options error: the memory expansion factor, expand, should be a single number!");
            opts.memoryExpansionFactor = (Real)val;
        }
        else if (!strcmp(fieldName, "reorth")) {
            if (nrows!=1 || ncols!=1 || val<0.0)
                mexErrMsgTxt("filtlan options error: the reorthogonalization method index, reorth, should be a single nonnegative number!");
            opts.reorth = (int)val;
        }
        else if (!strcmp(fieldName, "gamma1")) {
            if (nrows!=1 || ncols!=1)
                mexErrMsgTxt("filtlan options error: the local reorthogonalization coefficient, gamma1, should be a single number!");
            opts.localReorthGamma = (Real)val;
        }
        else if (!strcmp(fieldName, "gamma2")) {
            if (nrows!=1 || ncols!=1)
                mexErrMsgTxt("filtlan options error: the double reorthogonalization coefficient, gamma2, should be a single number!");
            opts.doubleReorthGamma = (Real)val;
        }
        else if (!strcmp(fieldName, "strategy")) {
            if (nrows!=1 || ncols!=1 || val<0.0)
                mexErrMsgTxt("filtlan options error: the partial reorthogonalization strategy index, strategy, should be a single nonnegative number!");
            opts.partialReorthStrategy = (int)val;
        }
        else if (!strcmp(fieldName, "checkorth")) {
            if (nrows!=1 || ncols!=1)
                mexErrMsgTxt("filtlan options error: the flag of whether to check semi-orthogonality or not, checkorth, should be a 0-1 boolean number!");
            opts.checkReorth = (bool)val;
        }
        else if (!strcmp(fieldName, "nev")) {
            if (nrows!=1 || ncols!=1 || val<1.0)
                mexErrMsgTxt("filtlan options error: the number of eigenvalues requested, nev, should be a single positive integer!");
            opts.neigWanted = (mkIndex)val;
        }
        else {
            mexPrintf("filtlan: invalid option '%s'!\n", fieldName);
            mexErrMsgTxt("There exists an invalid option!");
        }
    }
}

mxArray *convertFilteredLanczosInfo(const FilteredLanczosInfo &info) {
    // create a mxArray as a structure array of numFields fields
    const int numFields = 14;  // number of fields
    static const char fieldNames[numFields][32] = {
        // number of Lanczos iterations
        "numIter",
        // CPU time
        "forEigenRangeCpuTime",
        "forNextKrylovVectorCpuTime",
        "reorthogonalizationCpuTime",
        "convergenceCheckCpuTime",
        "forEigenVectorsCpuTime",
        "LanczosCpuTime",
        // memory required for the Lanczos iterations
        "memoryForLanczosInBytes",
        // whether all the eigenvalues are checked converged
        "allEigenvaluesCheckedConverged",
        // reorthogonalization cost
        "reorthIterCount",
        "reorthVectorCount",
        "reorthVectorRate",
        "localReorthIterCount",
        "doubleReorthIterCount",
    };
    // the corresponding class IDs
    static const mxClassID fieldClassIDs[numFields] = {
        mxUINT32_CLASS,
        #ifdef USE_SINGLE
            mxSINGLE_CLASS,
            mxSINGLE_CLASS,
            mxSINGLE_CLASS,
            mxSINGLE_CLASS,
            mxSINGLE_CLASS,
            mxSINGLE_CLASS,
        #else
            mxDOUBLE_CLASS,
            mxDOUBLE_CLASS,
            mxDOUBLE_CLASS,
            mxDOUBLE_CLASS,
            mxDOUBLE_CLASS,
            mxDOUBLE_CLASS,
        #endif
        mxUINT32_CLASS,
        mxUINT8_CLASS,
        mxUINT32_CLASS,
        mxUINT32_CLASS,
        #ifdef USE_SINGLE
            mxSINGLE_CLASS,
        #else
            mxDOUBLE_CLASS,
        #endif
        mxUINT32_CLASS,
        mxUINT32_CLASS
    };
    const char* fnames[numFields];
    for (int i=0; i<numFields; i++)
        fnames[i] = &(fieldNames[i][0]);
    mxArray *mxOut = mxCreateStructMatrix(1, 1, numFields, fnames);

    const void* ptr[numFields];
    ptr[0]  = &info.numIter;
    ptr[1]  = &info.forEigenRangeCpuTime;
    ptr[2]  = &info.forNextKrylovVectorCpuTime;
    ptr[3]  = &info.reorthogonalizationCpuTime;
    ptr[4]  = &info.convergenceCheckCpuTime;
    ptr[5]  = &info.forEigenVectorsCpuTime;
    ptr[6]  = &info.LanczosCpuTime;

    ptr[7]  = &info.memoryForLanczosInBytes;
    ptr[8]  = &info.allEigenvaluesCheckedConverged;

    ptr[9]  = &info.reorthIterCount;
    ptr[10] = &info.reorthVectorCount;
    ptr[11] = &info.reorthVectorRate;
    ptr[12] = &info.localReorthIterCount;
    ptr[13] = &info.doubleReorthIterCount;
    // &info.intervals can also be set

    // now copy the data from input to output
    for (int i=0; i<numFields; i++) {
        // allocate a space for output
        const mwSize one = 1;
        mxArray *aout = mxCreateNumericArray(1, &one, fieldClassIDs[i], mxREAL);
        void *dout = mxGetData(aout);
        int sz = mxGetElementSize(aout);  // number of bytes to be copied

        // now copy the data
        memcpy(dout, ptr[i], sz);

        // set the field in output structure
        mxSetFieldByNumber(mxOut, 0, i, aout);
    }
    return mxOut;
}

void printUsage() {
    FilteredLanczosOptions opts0;  // for the default filtered Lanczos options
    mexPrintf("FILTLAN finds the eigenvalues and eigenvectors of a real symmetric matrix\n");
    mexPrintf("in a specified interval using a polynomial filtered Lanczos procedure with\n");
    mexPrintf("partial or full reorthogonalization.\n");
    mexPrintf("\n");
    mexPrintf("D = FILTLAN(A,[a,b],polydeg,basedeg) returns a vector D of the eigenvalues\n");
    mexPrintf("of A in [a,b], where A is symmetric and sparse.\n");
    mexPrintf("\n");
    mexPrintf("[V,D] = FILTLAN(A,[a,b],polydeg,basedeg) returns a diagonal matrix D with\n");
    mexPrintf("diagonal elements the eigenvalues of A in [a,b], and a matrix V whose\n");
    mexPrintf("columns are the corresponding eigenvectors.\n");
    mexPrintf("\n");
    mexPrintf("[V,D,info] = FILTLAN(A,[a,b],polydeg,basedeg) in addition returns some\n");
    mexPrintf("information:\n");
    mexPrintf("    info.numIter: number of Lanczos iterations\n");
    mexPrintf("    info.forEigenRangeCpuTime: CPU time for determining the eigen-range\n");
    mexPrintf("    info.forNextKrylovVectorCpuTime: CPU time for matrix-vector products\n");
    mexPrintf("    info.reorthogonalizationCpuTime: CPU time for reorthogonalization\n");
    mexPrintf("    info.convergenceCheckCpuTime: CPU time for convergence check\n");
    mexPrintf("    info.forEigenVectorCpuTime: CPU time for obtaining the eigenvectors\n");
    mexPrintf("    info.LanczosCpuTime: this is the CPU time for about the entire procedure\n");
    mexPrintf("        (CPU time is measured in seconds in all cases)\n");
    mexPrintf("    info.memoryForLanczosInBytes: memory in bytes used for the Lanczos\n");
    mexPrintf("        iterations\n");
    mexPrintf("    info.allEigenvaluesCheckedConverged: 0 if all eigenvalues converged;\n");
    mexPrintf("        otherwise not all converged\n");
    mexPrintf("    info.reorthIterCount: number of iterations in which the\n");
    mexPrintf("        reorthogonalization is performed\n");
    mexPrintf("    info.reorthVectorCount: number of Lanczos vectors which the\n");
    mexPrintf("        reorthogonalization is performed against to\n");
    mexPrintf("    info.reorthVectorRate: rate of Lanczos vectors reorthogonalization is\n");
    mexPrintf("        performed against to, compared to full reorthogonalization\n");
    mexPrintf("    info.localReorthIterCount: number of iterations in which the extended\n");
    mexPrintf("        local reorthogonalization is performed\n");
    mexPrintf("    info.doubleReorthIterCount: number of iterations in which the double\n");
    mexPrintf("        reorthogonalization is performed\n");
    mexPrintf("\n");
    mexPrintf("FILTLAN(A,[a,b],polydeg,basedeg,opts) specifies the options:\n");
    mexPrintf("    opts.rseed: random seed, a nonnegative integer\n");
    mexPrintf("        If the initial Lanczos vector is not provided, a random vector is\n");
    mexPrintf("        generated with the random seed.\n");
    mexPrintf("    opts.eigrangeiter: number of Lanczos iterations to determine an interval\n");
    mexPrintf("        which (tightly) contains the eigenvalues of A\n");
    mexPrintf("    opts.maxiter: maximum number of Lanczos iterations\n");
    mexPrintf("    opts.miniter: minimum number of Lanczos iterations\n");
    mexPrintf("    opts.extraiter: extra number of Lanczos iterations after deemed\n");
    mexPrintf("        converged, to avoid missing eigenvalues and to improve the accuracy\n");
    mexPrintf("    opts.stride: convergence is checked every opts.stride Lanczos iterations\n");
    mexPrintf("    opts.v0: starting Lanczos vector (default is randomly generated)\n");
    mexPrintf("    opts.tol: tolerance for convergence test (default is %e)\n", opts0.tol);
    mexPrintf("    opts.nev: number of eigenvalues to be sought");
    if (opts0.neigWanted == 0)
        mexPrintf("\n        (default is to find all eigenvalues in the specified interval)\n");
    else
        mexPrintf(" (default is %u)\n", opts0.neigWanted);
    mexPrintf("    opts.disp: diagnostic information display level, 0");
    if (opts0.disp == 0)
        mexPrintf(" (default)");
    mexPrintf(", 1");
    if (opts0.disp == 1)
        mexPrintf(" (default)");
    mexPrintf(", or 2");
    if (opts0.disp == 2)
        mexPrintf(" (default)");
    mexPrintf("\n");
    mexPrintf("        0 (no information), 1 (some information), or 2 (more information)\n");
    mexPrintf("    opts.expand: specifies how large the memory should be expanded when the\n");
    mexPrintf("        allocated memory is insufficient to store the Lanczos vectors\n");
    mexPrintf("        (default is %.2g)\n", opts0.memoryExpansionFactor);
    mexPrintf("    opts.reorth: 1 or 2 (default is %u)\n", opts0.reorth);
    mexPrintf("        1 for partial reorthogonalization\n");
    mexPrintf("        2 for full reorthogonalization\n");
    mexPrintf("    opts.gamma2: reorthogonalization is doubled in the iterations with\n");
    mexPrintf("        nrm < gamma2*nrm_old, where nrm_old and nrm are the norms of the\n");
    mexPrintf("        latest Lanczos vectors before and after reorthogonalization\n");
    mexPrintf("        -- gamma2<=0.0: reorthogonalization is never doubled\n");
    mexPrintf("        -- gamma2>=1.0: reorthogonalization is always doubled\n");
    mexPrintf("        (by default, gamma2==%f)\n", opts0.doubleReorthGamma);
    mexPrintf("\n");
    mexPrintf("  options for partial reorthogonalization:\n");
    mexPrintf("    opts.gamma1: extended local reorthogonalization is performed in the\n");
    mexPrintf("        iterations with beta[j-1] < gamma1*beta[j], where beta[j] is the\n");
    mexPrintf("        latest beta and beta[j-1] is the second to the latest beta\n");
    mexPrintf("        -- gamma1<=0.0: local reorthogonalization is never performed\n");
    mexPrintf("        -- gamma1>=1.0: local reorthogonalization is always performed\n");
    mexPrintf("        (by default, gamma1==%f)\n", opts0.localReorthGamma);
    mexPrintf("    opts.strategy: 0");
    if (opts0.partialReorthStrategy == 0)
        mexPrintf(" (default)");
    mexPrintf(", 1");
    if (opts0.partialReorthStrategy == 1)
        mexPrintf(" (default)");
    mexPrintf(", 2");
    if (opts0.partialReorthStrategy == 2)
        mexPrintf(" (default)");
    mexPrintf(", or 3");
    if (opts0.partialReorthStrategy == 3)
        mexPrintf(" (default)");
    mexPrintf("\n");
    mexPrintf("        0 for reorthogonalization against previous Lanczos vectors in\n");
    mexPrintf("            groups; each group is v(i),v(i+1),...,v(j) with omega all less\n");
    mexPrintf("            then eta, and there is v(k) < delta, i<=k<=j\n");
    mexPrintf("        1 for reorthogonalization against previous Lanczos vectors v(i)\n");
    mexPrintf("            with omega(i) < eta\n");
    mexPrintf("        2 for reorthogonalization against previous Lanczos vectors v(i),\n");
    mexPrintf("            v(i+1),...,v(j), such that i is the smallest index for\n");
    mexPrintf("            omega(i) < eta, and j is the largest index for omega(j) < eta\n");
    mexPrintf("        3 for reorthogonalization against all previous Lanczos vectors\n");
    mexPrintf("        Here omega(i) is the estimated orthogonal error of v(i) against the\n");
    mexPrintf("        previous Lanczos vectors. By default, eta=eps^0.75, delta=eps^0.5.\n");
    mexPrintf("    opts.checkreorth: whether to check whether semi-orthogonalization is\n");
    mexPrintf("        preserved with partial reorthogonalization or not (default is %d)\n", opts0.checkReorth);
    mexPrintf("        opts.checkreorth should not be set as 1 for large problems, since\n");
    mexPrintf("        it is expensive.\n");
    mexPrintf("\n");
}
