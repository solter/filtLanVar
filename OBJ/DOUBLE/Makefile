ROOT_DIR      = ../..
LIB_DIR       = $(ROOT_DIR)/LIB

include $(ROOT_DIR)/Make.inc

MATKIT_DIR    = $(ROOT_DIR)/MATKIT
MATKIT_INC    = $(MATKIT_DIR)/INCLUDE
MATKIT_HEADS  = $(MATKIT_INC)/matkitdef.h $(MATKIT_INC)/matkitfunc.h $(MATKIT_INC)/basics.h $(MATKIT_INC)/vector.h $(MATKIT_INC)/matrix.h $(MATKIT_INC)/symmatrix.h $(MATKIT_INC)/spmatrix.h

FILTLAN_DIR   = $(ROOT_DIR)
FILTLAN_INC   = $(FILTLAN_DIR)/INCLUDE
FILTLAN_SRC   = $(FILTLAN_DIR)/SRC


all: $(LIB_DIR)/libdfiltlan.a

# by default, we use double precision (not single precision) floating point numbers
# XFLAGS being nothing means using the default
XFLAGS =

# symmetric tridiagonal eigensolver and Lanczos eigensolver
symtrieig.o:  $(FILTLAN_SRC)/symtrieig.cpp $(FILTLAN_INC)/symtrieig.h  $(MATKIT_HEADS)
	$(CPP) $(CFLAGS) $(XFLAGS) -I$(MATKIT_INC) -I$(FILTLAN_INC) -c $(FILTLAN_SRC)/symtrieig.cpp

laneig.o:     $(FILTLAN_SRC)/laneig.cpp $(FILTLAN_INC)/laneig.h  $(FILTLAN_INC)/symtrieig.h  $(MATKIT_HEADS)
	$(CPP) $(CFLAGS) $(XFLAGS) -I$(MATKIT_INC) -I$(FILTLAN_INC) -c $(FILTLAN_SRC)/laneig.cpp


# filtered Lanczos eigensolver
polyfilt.o:   $(FILTLAN_SRC)/polyfilt.cpp $(FILTLAN_INC)/polyfilt.h  $(MATKIT_HEADS)
	$(CPP) $(CFLAGS) $(XFLAGS) -I$(MATKIT_INC) -I$(FILTLAN_INC) -c $(FILTLAN_SRC)/polyfilt.cpp

filtlan.o:    $(FILTLAN_SRC)/filtlan.cpp $(FILTLAN_INC)/filtlan.h  $(FILTLAN_INC)/polyfilt.h $(FILTLAN_INC)/laneig.h $(FILTLAN_INC)/symtrieig.h  $(MATKIT_HEADS)
	$(CPP) $(CFLAGS) $(XFLAGS) -I$(MATKIT_INC) -I$(FILTLAN_INC) -c $(FILTLAN_SRC)/filtlan.cpp


# the archieve
OBJS = symtrieig.o laneig.o polyfilt.o filtlan.o

$(LIB_DIR)/libdfiltlan.a: $(OBJS)
	if [ ! -d $(LIB_DIR) ]; then \
    mkdir $(LIB_DIR); \
fi
	$(AR) $(ARFLAGS) $(LIB_DIR)/libdfiltlan.a $(OBJS)


clean:
	rm -f *.o $(LIB_DIR)/libdfiltlan.a
