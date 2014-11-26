ROOT_DIR      = ..

include $(ROOT_DIR)/Make.inc
include Make.inc

LIB_DIR       = $(ROOT_DIR)/LIB

MATKIT_DIR    = $(ROOT_DIR)/MATKIT
MATKIT_INC    = $(MATKIT_DIR)/INCLUDE
MATKIT_HEADS  = $(MATKIT_INC)/matkitdef.h $(MATKIT_INC)/matkitfunc.h  $(MATKIT_INC)/basics.h $(MATKIT_INC)/vector.h $(MATKIT_INC)/matrix.h $(MATKIT_INC)/symmatrix.h $(MATKIT_INC)/spmatrix.h
MATKIT_LIB    = $(LIB_DIR)/libdmatkit.a

FILTLAN_DIR   = ..
FILTLAN_INC   = $(FILTLAN_DIR)/INCLUDE
LANCZOS_HEADS = $(FILTLAN_INC)/laneig.h $(FILTLAN_INC)/symtrieig.h
FILTLAN_LIB   = $(LIB_DIR)/libdfiltlan.a

# by default, we use double precision (not single precision) floating numbers
# XFLAGS being nothing means using the default
XFLAGS =

all: arpack_driver silaneig_driver


# ARPACK implementation of implicited restarted Lanczos w/ full reorthogonalization
arpack.o:          arpack.cpp arpack.h $(MATKIT_HEADS)
	$(CPP) $(CFLAGS) $(XFLAGS) -I$(MATKIT_INC) -c arpack.cpp

arpack_driver.o:   arpack_driver.cpp arpack.h $(MATKIT_HEADS)
	$(CPP) $(CFLAGS) $(XFLAGS) -I$(MATKIT_INC) -c arpack_driver.cpp

arpack_driver:     arpack_driver.o arpack.o $(MATKIT_LIB) $(ARPACK_LIB)
	$(CPP) -o arpack_driver arpack_driver.o arpack.o $(MATKIT_LIB) $(ARPACK_LIB) $(LOADFLAGS)


# shift-and-invert Lanczos w/ partial reorthogonalization, requiring UMFPACK for sparse LU
silaneig.o:        silaneig.cpp silaneig.h $(MATKIT_HEADS) $(LANCZOS_HEADS)
	$(CPP) $(CFLAGS) $(XFLAGS) $(UMFPACK_INCX) -I$(MATKIT_INC) -I$(FILTLAN_INC) -c silaneig.cpp

silaneig_driver.o: silaneig_driver.cpp silaneig.h $(MATKIT_HEADS) $(LANCZOS_HEADS)
	$(CPP) $(CFLAGS) $(XFLAGS) $(UMFPACK_INCX) -I$(MATKIT_INC) -I$(FILTLAN_INC) -c silaneig_driver.cpp

silaneig_driver:   silaneig_driver.o silaneig.o $(FILTLAN_LIB) $(MATKIT_LIB) $(UMFPACK_LIBX)
	$(CPP) -o silaneig_driver silaneig_driver.o silaneig.o $(FILTLAN_LIB) $(MATKIT_LIB) $(UMFPACK_LIBX) $(LOADFLAGS)


clean:
	rm -f *.o arpack_driver silaneig_driver
