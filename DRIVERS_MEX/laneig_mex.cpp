//============================================================================
// The laneig mex driver for OCTAVE/MATLAB users
//
// Reference:
// "A Filtered Lanczos Procedure for Extreme and Interior Eigenvalue Problems"
// H.-r. Fang and Y. Saad, University of Minnesota Technical Report, 2011.
//============================================================================

#include <stdio.h>
#include <string.h>
#include "matkit.h"
#include "laneig.h"
#include "mex.h"

void printUsage();
void parseLanczosOptions(LanczosOptions &opts, const mxArray *mxOpts);
mxArray *convertLanczosInfo(const LanczosInfo &info);


// OCTAVE/MATLAB: [V, D, info] = laneig(A, nev, part, opts);
//                [V, D]       = laneig(A, nev, part, opts);
//                    D        = laneig(A, nev, part, opts);
// ======
// input:
// A is a sparse matrix whose eigenvalues are of interest
// nev is the number of eigenvalues requested (default 6)
// part indicates which part of eigenvalues are sought:
//    "SA" - Smallest Algebraic, for the smallest eigenvalues
//    "LA" - Largest Algebraic, for the largest eigenvalues (default)
//    "BE" - Both Ends, one more from high end if nev is odd
// opts is a collection of Lanczos options, which is optional, and
// by default it is from the default constructor of LanczosOptions
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
    if (nrhs > 4)
        mexErrMsgTxt("laneig: at least 1 and at most 4 input arguments should be provided!");
    if (nlhs > 3)
        mexErrMsgTxt("laneig: at most 3 output arguments are generated!");

    // OCTAVE/MATLAB: [nrows, ncols] = size(A);
    mkIndex nrows = (mkIndex)mxGetM(prhs[0]);
    mkIndex ncols = (mkIndex)mxGetN(prhs[0]);
    if (nrows != ncols)
        mexErrMsgTxt("laneig: input matrix must be square!");

    if (!mxIsSparse(prhs[0]))
        mexErrMsgTxt("laneig: input as a dense matrix is not (yet) supported!");

    SparseMatrix A = SparseMatrix::mxArray2SparseMatrix(prhs[0]);
    // the routine mxArray2SparseMatrix is declared in spmatrix.h of MATKIT
    // assume that the input matrix A is symmetric

    // get the requested number of eigenvalues from the 2nd input argument, if any
    mkIndex neigWanted = 6;
    if (nrhs >= 2) {
        if (!mxIsNumeric(prhs[1]) || mxGetM(prhs[1])!=1 || mxGetN(prhs[1])!=1)
            mexErrMsgTxt("laneig: the 2nd input argument, the number of eigenvalues requested, should be a single positive integer!");
        double nev = mxArray_to_double(prhs[1]);
        if (nev < 1.0)
            mexErrMsgTxt("laneig: the 2nd input argument, the number of eigenvalues requested, should be a single positive integer!");
        neigWanted = (mkIndex)nev;
    }

    // determine which part of eigenvalues requested from the 3rd input argument, if any
    char defaultEigPart[] = "LA";
    char *eigPart = defaultEigPart;
    if (nrhs >= 3) {
        if (!mxIsChar(prhs[2]))
            mexErrMsgTxt("laneig: the 3rd input argument should consist of two characters, 'LA', 'SA', or 'BE', specifying the part of eigenvalues requested!");
        eigPart = mxArrayToString(prhs[2]);
    }

    // set options from the 4th input argument, if any
    LanczosOptions opts;  // the default is set by the constructor
    if (nrhs == 4) {
        if (!mxIsStruct(prhs[3])) {
            mexErrMsgTxt("laneig: the 4th input argument must be structure array for Lanczos options!");
        }
        parseLanczosOptions(opts, prhs[3]);
    }

    if (nlhs <= 1)
        opts.wantEigVec = false;  // the 1st output argument is a vector containing eigenvalues; eigenvectors are not requested
    else
        opts.wantEigVec = true;   // the 1st output argument is a matrix with columns as eigenvectors; the 2nd output argument is a diagonal sparse matrix of eigenvalues

    // find the eigenvalues/eigenvectors requested
    Vector lambda;
    Matrix V;
    LanczosInfo outputInfo = LanczosEigenSolver(lambda, V, A, neigWanted, eigPart, opts);
    if (nrhs >= 3)
        mxFree(eigPart);  // the output of mxArrayToString should be mxFree-ed

    // set the output
    if (nlhs <= 1) {
        plhs[0] = Vector::Vector2mxArray(lambda);
    }
    else {  // nlhs >= 2
        plhs[0] = Matrix::Matrix2mxArray(V);  // 1st output argument
        plhs[1] = SparseMatrix::SparseMatrix2mxArray(lambda.spdiag());  // 2nd output argument
        if (nlhs == 3) {
            // now set the 3rd output argument
            plhs[2] = convertLanczosInfo(outputInfo);
        }
    }
}

void parseLanczosOptions(LanczosOptions &opts, const mxArray *mxOpts) {
    // get input arguments
    int numFields = mxGetNumberOfFields(mxOpts);
    if (mxGetNumberOfElements(mxOpts) != 1) {
        mexErrMsgTxt("laneig options parser: each field should have one and only one element!");
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
        double val = mxArray_to_double(fieldData);
        if (!strcmp(fieldName, "rseed")) {
            if (nrows!=1 || ncols!=1 || val<0.0)
                mexErrMsgTxt("laneig options error: the random seed, rseed, should be a single nonnegative integer!");
            srand((unsigned)val);
        }
        else if (!strcmp(fieldName, "maxiter")) {
            if (nrows!=1 || ncols!=1 || val<0.0)
                mexErrMsgTxt("laneig options error: maximum number of iterations, maxiter, should be a single nonnegative number!");
            opts.maxIter = (mkIndex)val;
        }
        else if (!strcmp(fieldName, "miniter")) {
            if (nrows!=1 || ncols!=1 || val<0.0)
                mexErrMsgTxt("laneig options error: minimum number of iterations, miniter, should be a single nonnegative number!");
            opts.minIter = (mkIndex)val;
        }
        else if (!strcmp(fieldName, "extraiter")) {
            if (nrows!=1 || ncols!=1 || val<0.0)
                mexErrMsgTxt("laneig options error: extra number of iterations (after deemed convergence), extraiter, should be a single nonnegative number!");
            opts.extraIter = (mkIndex)val;
        }
        else if (!strcmp(fieldName, "stride")) {
            if (nrows!=1 || ncols!=1 || val<0.0)
                mexErrMsgTxt("laneig options error: the convergence to be checked every 'stride' iterations, the stride should be a single nonnegative number!");
            opts.stride = (mkIndex)val;
        }
        else if (!strcmp(fieldName, "tol")) {
            if (nrows!=1 || ncols!=1)
                mexErrMsgTxt("laneig options error: the tolerance for checking convergence of eigenvalues, tol, should be a single number!");
            opts.tol = (Real)val;
        }
        else if (!strcmp(fieldName, "cut0")) {
            if (nrows!=1 || ncols!=1)
                mexErrMsgTxt("laneig options error: the upper bound of the eigenvalues to be sought in the lower end of the spectrum, cut0, should be a single number!");
            opts.eigLowCut = (Real)val;
        }
        else if (!strcmp(fieldName, "cut1")) {
            if (nrows!=1 || ncols!=1)
                mexErrMsgTxt("laneig options error: the lower bound of the eigenvalues to be sought in the upper end of the spectrum, cut1, should be a single number!");
            opts.eigHighCut = (Real)val;
        }
        else if (!strcmp(fieldName, "cut")) {
            // it means to set eigLowCut = eigHighCut = cut
            if (nrows!=1 || ncols!=1)
                mexErrMsgTxt("laneig options error: the bound of the eigenvalues, cut, should be a single number!");
            opts.eigLowCut = opts.eigHighCut = (Real)val;
        }
        else if (!strcmp(fieldName, "disp")) {
            if (nrows!=1 || ncols!=1 || val<0.0)
                mexErrMsgTxt("laneig options error: the diagnostic information display level, disp, should be a single nonnegative number!");
            opts.disp = (int)val;
        }
        else if (!strcmp(fieldName, "expand")) {
            if (nrows!=1 || ncols!=1)
                mexErrMsgTxt("laneig options error: the memory expansion factor, expand, should be a single number!");
            opts.memoryExpansionFactor = (Real)val;
        }
        else if (!strcmp(fieldName, "reorth")) {
            if (nrows!=1 || ncols!=1 || val<0.0)
                mexErrMsgTxt("laneig options error: the reorthogonalization method index, reorth, should be a single nonnegative number!");
            opts.reorth = (int)val;
        }
        else if (!strcmp(fieldName, "gamma1")) {
            if (nrows!=1 || ncols!=1)
                mexErrMsgTxt("laneig options error: the local reorthogonalization coefficient, gamma1, should be a single number!");
            opts.localReorthGamma = (Real)val;
        }
        else if (!strcmp(fieldName, "gamma2")) {
            if (nrows!=1 || ncols!=1)
                mexErrMsgTxt("laneig options error: the double reorthogonalization coefficient, gamma2, should be a single number!");
            opts.doubleReorthGamma = (Real)val;
        }
        else if (!strcmp(fieldName, "strategy")) {
            if (nrows!=1 || ncols!=1 || val<0.0)
                mexErrMsgTxt("laneig options error: the partial reorthogonalization strategy index, strategy, should be a single nonnegative number!");
            opts.partialReorthStrategy = (int)val;
        }
        else if (!strcmp(fieldName, "checkorth")) {
            if (nrows!=1 || ncols!=1)
                mexErrMsgTxt("laneig options error: the flag of whether to check semi-orthogonality or not, checkorth, should be a 0-1 boolean number!");
            opts.checkReorth = (bool)val;
        }
        else {
            mexPrintf("laneig: invalid option '%s'!\n", fieldName);
            mexErrMsgTxt("There exists an invalid option!");
        }
    }
}

mxArray *convertLanczosInfo(const LanczosInfo &info) {
    // create a mxArray as a structure array of numFields fields
    const int numFields = 13;  // number of fields
    static const char fieldNames[numFields][32] = {
        // number of Lanczos iterations
        "numIter",
        // CPU time
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
        #else
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
    ptr[1]  = &info.forNextKrylovVectorCpuTime;
    ptr[2]  = &info.reorthogonalizationCpuTime;
    ptr[3]  = &info.convergenceCheckCpuTime;
    ptr[4]  = &info.forEigenVectorsCpuTime;
    ptr[5]  = &info.LanczosCpuTime;

    ptr[6]  = &info.memoryForLanczosInBytes;
    ptr[7]  = &info.allEigenvaluesCheckedConverged;

    ptr[8]  = &info.reorthIterCount;
    ptr[9]  = &info.reorthVectorCount;
    ptr[10] = &info.reorthVectorRate;
    ptr[11] = &info.localReorthIterCount;
    ptr[12] = &info.doubleReorthIterCount;

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
    LanczosOptions opts0;  // for the default Lanczos options
    mexPrintf("LANEIG finds a small number of eigenvalues and eigenvectors of a real\n");
    mexPrintf("symmetric matrix using the Lanczos method with partial or full\n");
    mexPrintf("reorthogonalization and without restarting.\n");
    mexPrintf("\n");
    mexPrintf("D = LANEIG(A) returns a vector of the 6 largest eigenvalues of A, where A\n");
    mexPrintf("is symmetric and sparse.\n");
    mexPrintf("\n");
    mexPrintf("[V,D] = LANEIG(A) returns a diagonal matrix D of the 6 largest eigenvalues\n");
    mexPrintf("of A, and a matrix V whose columns are the corresponding eigenvectors.\n");
    mexPrintf("\n");
    mexPrintf("[V,D,info] = LANEIG(A) in addition returns some information of the process:\n");
    mexPrintf("    info.numIter: number of Lanczos iterations\n");
    mexPrintf("    info.forNextKrylovVectorCpuTime: CPU time for matrix-vector products\n");
    mexPrintf("    info.reorthogonalizationCpuTime: CPU time for reorthogonalization\n");
    mexPrintf("    info.convergenceCheckCpuTime: CPU time for convergence check\n");
    mexPrintf("    info.forEigenVectorCpuTime: CPU time for obtaining the eigenvectors\n");
    mexPrintf("    info.LanczosCpuTime: this is the CPU time for about the entire procedure\n");
    mexPrintf("        (CPU time is measured in seconds in all cases)\n");
    mexPrintf("    info.memoryForLanczosInBytes: memory in bytes required for the Lanczos\n");
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
    mexPrintf("LANEIG(A,nev) returns the nev largest eigenvalues.\n");
    mexPrintf("\n");
    mexPrintf("LANEIG(A,nev,part) returns nev eigenvalues; part may be 'LA', 'SA', or 'BE':\n");
    mexPrintf("    'LA' for Largest Algebraic\n");
    mexPrintf("    'SA' for Smallest Algebraic\n");
    mexPrintf("    'BE' for Both Ends, one more from high end if nev is odd\n");
    mexPrintf("\n");
    mexPrintf("LANEIG(A,nev,part,opts) specifies the options:\n");
    mexPrintf("    opts.rseed: random seed, a nonnegative integer\n");
    mexPrintf("        If the initial Lanczos vector opts.v0 is not provided, a random\n");
    mexPrintf("        vector is generated with the random seed.\n");
    mexPrintf("    opts.maxiter: maximum number of Lanczos iterations\n");
    mexPrintf("    opts.miniter: minimum number of Lanczos iterations\n");
    mexPrintf("    opts.extraiter: extra number of Lanczos iterations after deemed\n");
    mexPrintf("        converged, to improve the accuracy and to avoid missing eigenvalues\n");
    mexPrintf("    opts.stride: convergence is checked every opts.stride Lanczos iterations\n");
    mexPrintf("    opts.v0: starting Lanczos vector (default is randomly generated)\n");
    mexPrintf("    opts.tol: tolerance for convergence test (default is %e)\n", opts0.tol);
    mexPrintf("    opts.cut0: upper bound of the eigenvalues to be sought in the lower end\n");
    mexPrintf("        of the spectrum; effective if the part of eigenvalues requested is\n");
    mexPrintf("        'SA' or 'BE'; default is no bound\n");
    mexPrintf("    opts.cut1: lower bound of the eigenvalues to be sought in the upper end\n");
    mexPrintf("        of the spectrum; effective if the part of eigenvalues requested is\n");
    mexPrintf("        'LA' or 'BE'; default is no bound\n");
    mexPrintf("    opts.cut: set both opts.cut0 and opts.cut1 as the same value opts.cut\n");
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
    mexPrintf("        0 for no information; 1 for some information; 2 for more information\n");
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
    mexPrintf("        (by default, gamma2==%g)\n", opts0.doubleReorthGamma);
    mexPrintf("\n");
    mexPrintf("  options for partial reorthogonalization:\n");
    mexPrintf("    opts.gamma1: extended local reorthogonalization is performed in the\n");
    mexPrintf("        iterations with beta[j-1] < gamma1*beta[j], where beta[j] is the\n");
    mexPrintf("        latest beta and beta[j-1] is the second to the latest beta\n");
    mexPrintf("        -- gamma1<=0.0: local reorthogonalization is never performed\n");
    mexPrintf("        -- gamma1>=1.0: local reorthogonalization is always performed\n");
    mexPrintf("        (by default, gamma1==%g)\n", opts0.localReorthGamma);
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
