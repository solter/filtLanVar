#-----------------------------------------------------------------------------
# compilers
#-----------------------------------------------------------------------------

CC = gcc
CPP = g++

# for OCTAVE mkoctave compilation, use "MEX_FLAG = FOR_OCTAVE"
# for MATLAB mex compilation, use "MEX_FLAG = FOR_MATLAB"
# for no mex compliation, use "MEX_FLAG ="
MEX_FLAG = FOR_OCTAVE
#MEX_FLAG = FOR_MATLAB
#MEX_FLAG =

ifeq ($(MEX_FLAG),FOR_MATLAB)
    # for MATLAB
    MATLABROOT = /soft/matlab2011b
    MEX = $(MATLABROOT)/bin/mex
else
    # for OCTAVE
    MEX = mkoctfile
endif


#-----------------------------------------------------------------------------
# flags for compilers
#-----------------------------------------------------------------------------

#CFLAGS = -g -O0 -Wall -DUSE_NAMESPACE -DUSE_BLAS
CFLAGS = -g -O0 -Wall -DUSE_NAMESPACE
LOADFLAGS = -llapack -lblas
## the following line is for "link" in itasca.msi.umn.edu, assuming "-DUSE_BLAS" and/or "-DUSE_SPBLAS" is set in $(CFLAGS)
# sequential
#LOADFLAGS = -L/soft/intel/mkl/10.2.1.017/lib/em64t -lmkl_lapack -lmkl_intel_lp64 -lmkl_sequential -lmkl_core
# multi-core threading
#LOADFLAGS = -L/soft/intel/mkl/10.2.1.017/lib/em64t -lmkl_lapack -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lguide -lpthread -lm
#LOADFLAGS = -L/soft/intel/mkl/10.2.1.017/lib/em64t -lmkl_lapack -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core -liomp5 -lpthread -lm

ifeq ($(MEX_FLAG),FOR_OCTAVE)
    # flags for the OCTAVE compiler mkoctfile
    MEXCFLAGS =    -g -DUSE_MEX -DUSE_BLAS
    # MEXCFLAGS =    -DUSE_MEX -DUSE_BLAS
    MEXLOADFLAGS = -llapack -lblas
    # MEXEXT is the extension of compiled mex-files
    # it is also the output "mexext" in OCTAVE console
    MEXEXT = mex
else ifeq ($(MEX_FLAG),FOR_MATLAB)
    # flags for the MATLAB compiler mex
    # "-largeArrayDims" is required for 64-bit OS machines
    #MEXCFLAGS =    -cxx -g -largeArrayDims -DUSE_MEX -DUSE_MWBLAS -DUSE_MWLAPACK
MEXCFLAGS =    -cxx -g -largeArrayDims -DUSE_MEX -DUSE_MWLAPACK
    # MEXCFLAGS =    -cxx -O -largeArrayDims -DUSE_MEX -DUSE_MWBLAS -DUSE_MWLAPACK
    MEXLOADFLAGS = -cxx -llapack -lblas
    # MEXEXT is the extension of compiled mex-files
    # it is also the output of $(MATLABROOT)/bin/mexext, and also the output of "mex" in MATLAB console
    MEXEXT = mexa64
endif


#-----------------------------------------------------------------------------
# archive command
#-----------------------------------------------------------------------------

AR = ar
ARFLAGS = cr
