ROOT_DIR      = ..
LIB_DIR       = $(ROOT_DIR)/LIB

include $(ROOT_DIR)/Make.inc

MATKIT_DIR    = $(ROOT_DIR)/MATKIT
MATKIT_INC    = $(MATKIT_DIR)/INCLUDE
MATKIT_HEADS  = $(MATKIT_INC)/matkitdef.h $(MATKIT_INC)/matkitfunc.h  $(MATKIT_INC)/basics.h $(MATKIT_INC)/vector.h $(MATKIT_INC)/matrix.h $(MATKIT_INC)/symmatrix.h $(MATKIT_INC)/spmatrix.h
MATKIT_LIB    = $(LIB_DIR)/libdmatkit_mex.a

FILTLAN_DIR   = $(ROOT_DIR)
FILTLAN_INC   = $(FILTLAN_DIR)/INCLUDE
FILTLAN_HEADS = $(FILTLAN_INC)/filtlan.h $(FILTLAN_INC)/polyfilt.h $(FILTLAN_INC)/laneig.h $(FILTLAN_INC)/symtrieig.h
FILTLAN_LIB   = $(LIB_DIR)/libdfiltlan_mex.a

# by default, we use double precision (not single precision) floating point numbers
# XFLAGS being nothing means using the default
XFLAGS =

DRIVERS_SRC = .

DRIVERS = laneig.$(MEXEXT) filtlan.$(MEXEXT)

ifeq ($(MEX_FLAG),FOR_OCTAVE)
    all: $(DRIVERS)
else ifeq ($(MEX_FLAG),FOR_MATLAB)
    all: $(DRIVERS)
else
    all:
endif

ifeq ($(MEX_FLAG),FOR_OCTAVE)
# for OCTAVE mex files
laneig.$(MEXEXT):   laneig.o $(FILTLAN_LIB) $(MATKIT_LIB)
	$(MEX) laneig.o $(FILTLAN_LIB) $(MATKIT_LIB) $(MEXLOADFLAGS)
	mv laneig.oct laneig.$(MEXEXT)

laneig.o:       $(DRIVERS_SRC)/laneig_mex.cpp $(MATKIT_HEADS) $(FILTLAN_INC)/laneig.h $(FILTLAN_INC)/symtrieig.h
	$(MEX) $(MEXCFLAGS) $(XFLAGS) -I$(MATKIT_INC) -I$(FILTLAN_INC) -c $(DRIVERS_SRC)/laneig_mex.cpp -o laneig_mex.o
	mv laneig_mex.o laneig.o

filtlan.$(MEXEXT):  filtlan.o $(FILTLAN_LIB) $(MATKIT_LIB)
	$(MEX) filtlan.o $(FILTLAN_LIB) $(MATKIT_LIB) $(MEXLOADFLAGS)
	mv filtlan.oct filtlan.$(MEXEXT)

filtlan.o:      $(DRIVERS_SRC)/filtlan_mex.cpp $(MATKIT_HEADS) $(FILTLAN_HEADS)
	$(MEX) $(MEXCFLAGS) $(XFLAGS) -I$(MATKIT_INC) -I$(FILTLAN_INC) -c $(DRIVERS_SRC)/filtlan_mex.cpp -o filtlan_mex.o
	mv filtlan_mex.o filtlan.o
else ifeq ($(MEX_FLAG),FOR_MATLAB)
# for MATLAB mex files
laneig.$(MEXEXT):   laneig_mex.o $(FILTLAN_LIB) $(MATKIT_LIB)
	$(MEX) -o laneig.$(MEXEXT) laneig_mex.o $(FILTLAN_LIB) $(MATKIT_LIB) $(MEXLOADFLAGS)

laneig_mex.o:       $(DRIVERS_SRC)/laneig_mex.cpp $(MATKIT_HEADS) $(FILTLAN_INC)/laneig.h $(FILTLAN_INC)/symtrieig.h
	$(MEX) $(MEXCFLAGS) $(XFLAGS) -I$(MATKIT_INC) -I$(FILTLAN_INC) -c $(DRIVERS_SRC)/laneig_mex.cpp

filtlan.$(MEXEXT):  filtlan_mex.o $(FILTLAN_LIB) $(MATKIT_LIB)
	$(MEX) -o filtlan.$(MEXEXT) filtlan_mex.o $(FILTLAN_LIB) $(MATKIT_LIB) $(MEXLOADFLAGS)

filtlan_mex.o:      $(DRIVERS_SRC)/filtlan_mex.cpp $(MATKIT_HEADS) $(FILTLAN_HEADS)
	$(MEX) $(MEXCFLAGS) $(XFLAGS) -I$(MATKIT_INC) -I$(FILTLAN_INC) -c $(DRIVERS_SRC)/filtlan_mex.cpp
endif


clean:
ifneq ($(MEX_FLAG),NO_MEX)
	rm -f *.o $(DRIVERS)
endif
