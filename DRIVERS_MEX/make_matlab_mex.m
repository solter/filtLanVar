% make the laneig mex file
disp('mex -cxx -largeArrayDims -DUSE_MEX -DUSE_MWBLAS -DUSE_MWLAPACK -I../MATKIT/INCLUDE -I../INCLUDE -o laneig laneig_mex.cpp ../SRC/laneig.cpp ../SRC/symtrieig.cpp ../MATKIT/SRC/matkitfunc.cpp ../MATKIT/SRC/basics.cpp ../MATKIT/SRC/vector.cpp ../MATKIT/SRC/matrix.cpp ../MATKIT/SRC/symmatrix.cpp ../MATKIT/SRC/spmatrix.cpp -lblas -llapack');
mex -cxx -largeArrayDims -DUSE_MEX -DUSE_MWBLAS -DUSE_MWLAPACK -I../MATKIT/INCLUDE -I../INCLUDE -o laneig laneig_mex.cpp ../SRC/laneig.cpp ../SRC/symtrieig.cpp ../MATKIT/SRC/matkitfunc.cpp ../MATKIT/SRC/basics.cpp ../MATKIT/SRC/vector.cpp ../MATKIT/SRC/matrix.cpp ../MATKIT/SRC/symmatrix.cpp ../MATKIT/SRC/spmatrix.cpp -lblas -llapack

% make the filtlan mex file
disp('mex -cxx -largeArrayDims -DUSE_MEX -DUSE_MWBLAS -DUSE_MWLAPACK -I../MATKIT/INCLUDE -I../INCLUDE -o filtlan filtlan_mex.cpp ../SRC/filtlan.cpp ../SRC/polyfilt.cpp ../SRC/laneig.cpp ../SRC/symtrieig.cpp ../MATKIT/SRC/matkitfunc.cpp ../MATKIT/SRC/basics.cpp ../MATKIT/SRC/vector.cpp ../MATKIT/SRC/matrix.cpp ../MATKIT/SRC/symmatrix.cpp ../MATKIT/SRC/spmatrix.cpp -lblas -llapack');
mex -cxx -largeArrayDims -DUSE_MEX -DUSE_MWBLAS -DUSE_MWLAPACK -I../MATKIT/INCLUDE -I../INCLUDE -o filtlan filtlan_mex.cpp ../SRC/filtlan.cpp ../SRC/polyfilt.cpp ../SRC/laneig.cpp ../SRC/symtrieig.cpp ../MATKIT/SRC/matkitfunc.cpp ../MATKIT/SRC/basics.cpp ../MATKIT/SRC/vector.cpp ../MATKIT/SRC/matrix.cpp ../MATKIT/SRC/symmatrix.cpp ../MATKIT/SRC/spmatrix.cpp -lblas -llapack

% remove the object files
disp('system(''rm -f *.o ../SRC/*.o ../MATKIT/SRC/*.o'');');
system('rm -f *.o ../SRC/*.o ../MATKIT/SRC/*.o');
