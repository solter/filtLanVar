ROOT_DIR      = ../..
LIB_DIR       = $(ROOT_DIR)/LIB

include $(ROOT_DIR)/Make.inc

MATKIT_DIR    = $(ROOT_DIR)/MATKIT
MATKIT_INC    = $(MATKIT_DIR)/INCLUDE
MATKIT_HEADS  = $(MATKIT_INC)/matkitdef.h $(MATKIT_INC)/matkitfunc.h $(MATKIT_INC)/basics.h $(MATKIT_INC)/vector.h $(MATKIT_INC)/matrix.h $(MATKIT_INC)/symmatrix.h $(MATKIT_INC)/spmatrix.h

FILTLAN_DIR   = $(ROOT_DIR)
FILTLAN_INC   = $(FILTLAN_DIR)/INCLUDE
FILTLAN_SRC   = $(FILTLAN_DIR)/SRC


ifeq ($(MEX_FLAG),FOR_OCTAVE)
all: $(LIB_DIR)/libsfiltlan_mex.a
else ifeq ($(MEX_FLAG),FOR_MATLAB)
all: $(LIB_DIR)/libsfiltlan_mex.a
else
all:
endif

XFLAGS = -DUSE_SINGLE

ifeq ($(MEX_FLAG),FOR_OCTAVE)
## for OCTAVE mex objective files
# symmetric tridiagonal eigensolver and Lanczos eigensolver
symtrieig.o:  $(FILTLAN_SRC)/symtrieig.cpp $(FILTLAN_INC)/symtrieig.h  $(MATKIT_HEADS)
	$(MEX) $(MEXCFLAGS) $(XFLAGS) -I$(MATKIT_INC) -I$(FILTLAN_INC) -c $(FILTLAN_SRC)/symtrieig.cpp -o symtrieig.o

laneig.o:     $(FILTLAN_SRC)/laneig.cpp $(FILTLAN_INC)/laneig.h  $(FILTLAN_INC)/symtrieig.h  $(MATKIT_HEADS)
	$(MEX) $(MEXCFLAGS) $(XFLAGS) -I$(MATKIT_INC) -I$(FILTLAN_INC) -c $(FILTLAN_SRC)/laneig.cpp -o laneig.o

# filtered Lanczos eigensolver
polyfilt.o:   $(FILTLAN_SRC)/polyfilt.cpp $(FILTLAN_INC)/polyfilt.h  $(MATKIT_HEADS)
	$(MEX) $(MEXCFLAGS) $(XFLAGS) -I$(MATKIT_INC) -I$(FILTLAN_INC) -c $(FILTLAN_SRC)/polyfilt.cpp -o polyfilt.o

filtlan.o:    $(FILTLAN_SRC)/filtlan.cpp $(FILTLAN_INC)/filtlan.h  $(FILTLAN_INC)/polyfilt.h $(FILTLAN_INC)/laneig.h $(FILTLAN_SRC)/symtrieig.cpp  $(MATKIT_HEADS)
	$(MEX) $(MEXCFLAGS) $(XFLAGS) -I$(MATKIT_INC) -I$(FILTLAN_INC) -c $(FILTLAN_SRC)/filtlan.cpp -o filtlan.o
else ifeq ($(MEX_FLAG),FOR_MATLAB)
## for MATLAB mex objective files
# symmetric tridiagonal eigensolver and Lanczos eigensolver
symtrieig.o:  $(FILTLAN_SRC)/symtrieig.cpp $(FILTLAN_INC)/symtrieig.h  $(MATKIT_HEADS)
	$(MEX) $(MEXCFLAGS) $(XFLAGS) -I$(MATKIT_INC) -I$(FILTLAN_INC) -c $(FILTLAN_SRC)/symtrieig.cpp

laneig.o:     $(FILTLAN_SRC)/laneig.cpp $(FILTLAN_INC)/laneig.h  $(FILTLAN_INC)/symtrieig.h  $(MATKIT_HEADS)
	$(MEX) $(MEXCFLAGS) $(XFLAGS) -I$(MATKIT_INC) -I$(FILTLAN_INC) -c $(FILTLAN_SRC)/laneig.cpp

# filtered Lanczos eigensolver
polyfilt.o:   $(FILTLAN_SRC)/polyfilt.cpp $(FILTLAN_INC)/polyfilt.h  $(MATKIT_HEADS)
	$(MEX) $(MEXCFLAGS) $(XFLAGS) -I$(MATKIT_INC) -I$(FILTLAN_INC) -c $(FILTLAN_SRC)/polyfilt.cpp

filtlan.o:    $(FILTLAN_SRC)/filtlan.cpp $(FILTLAN_INC)/filtlan.h  $(FILTLAN_INC)/polyfilt.h $(FILTLAN_INC)/laneig.h $(FILTLAN_SRC)/symtrieig.cpp  $(MATKIT_HEADS)
	$(MEX) $(MEXCFLAGS) $(XFLAGS) -I$(MATKIT_INC) -I$(FILTLAN_INC) -c $(FILTLAN_SRC)/filtlan.cpp
endif


# the archieve
OBJS = symtrieig.o laneig.o polyfilt.o filtlan.o

$(LIB_DIR)/libsfiltlan_mex.a: $(OBJS)
	if [ ! -d $(LIB_DIR) ]; then \
    mkdir $(LIB_DIR); \
fi
	$(AR) $(ARFLAGS) $(LIB_DIR)/libsfiltlan_mex.a $(OBJS)


clean:
ifneq ($(MEX_FLAG),NO_MEX)
	rm -f *.o $(LIB_DIR)/libsfiltlan_mex.a
endif
