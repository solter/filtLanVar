ROOT_DIR   = ../../..
MATKIT_DIR = ../..

include $(ROOT_DIR)/Make.inc

INC_DIR    = $(MATKIT_DIR)/INCLUDE
SRC_DIR    = $(MATKIT_DIR)/SRC
LIB_DIR    = $(ROOT_DIR)/LIB


ifeq ($(MEX_FLAG),FOR_OCTAVE)
all: $(LIB_DIR)/libsmatkit_mex.a
else ifeq ($(MEX_FLAG),FOR_MATLAB)
all: $(LIB_DIR)/libsmatkit_mex.a
else
all:
endif


XFLAGS = -DUSE_SINGLE


# the fundamental classes
HEADERS = $(INC_DIR)/matkitdef.h $(INC_DIR)/matkitfunc.h  $(INC_DIR)/basics.h $(INC_DIR)/vector.h $(INC_DIR)/matrix.h $(INC_DIR)/symmatrix.h $(INC_DIR)/spmatrix.h

ifeq ($(MEX_FLAG),FOR_OCTAVE)
# for OCTAVE mex object files
matkitfunc.o:  $(INC_DIR)/matkitdef.h $(INC_DIR)/matkitfunc.h $(SRC_DIR)/matkitfunc.cpp
	$(MEX) $(MEXCFLAGS) $(XFLAGS) -I$(INC_DIR) -c $(SRC_DIR)/matkitfunc.cpp -o matkitfunc.o

basics.o:      $(INC_DIR)/matkitdef.h $(INC_DIR)/basics.h $(SRC_DIR)/basics.cpp
	$(MEX) $(MEXCFLAGS) $(XFLAGS) -I$(INC_DIR) -c $(SRC_DIR)/basics.cpp -o basics.o

vector.o:      $(HEADERS) $(SRC_DIR)/vector.cpp
	$(MEX) $(MEXCFLAGS) $(XFLAGS) -I$(INC_DIR) -c $(SRC_DIR)/vector.cpp -o vector.o

matrix.o:      $(HEADERS) $(SRC_DIR)/matrix.cpp
	$(MEX) $(MEXCFLAGS) $(XFLAGS) -I$(INC_DIR) -c $(SRC_DIR)/matrix.cpp -o matrix.o

symmatrix.o:   $(HEADERS) $(SRC_DIR)/symmatrix.cpp
	$(MEX) $(MEXCFLAGS) $(XFLAGS) -I$(INC_DIR) -c $(SRC_DIR)/symmatrix.cpp -o symmatrix.o

spmatrix.o:    $(HEADERS) $(SRC_DIR)/spmatrix.cpp
	$(MEX) $(MEXCFLAGS) $(XFLAGS) -I$(INC_DIR) -c $(SRC_DIR)/spmatrix.cpp -o spmatrix.o
else ifeq ($(MEX_FLAG),FOR_MATLAB)
# for MATLAB mex object files
matkitfunc.o:  $(INC_DIR)/matkitdef.h $(INC_DIR)/matkitfunc.h $(SRC_DIR)/matkitfunc.cpp
	$(MEX) $(MEXCFLAGS) $(XFLAGS) -I$(INC_DIR) -c $(SRC_DIR)/matkitfunc.cpp

basics.o:      $(INC_DIR)/matkitdef.h $(INC_DIR)/basics.h $(SRC_DIR)/basics.cpp
	$(MEX) $(MEXCFLAGS) $(XFLAGS) -I$(INC_DIR) -c $(SRC_DIR)/basics.cpp

vector.o:      $(HEADERS) $(SRC_DIR)/vector.cpp
	$(MEX) $(MEXCFLAGS) $(XFLAGS) -I$(INC_DIR) -c $(SRC_DIR)/vector.cpp

matrix.o:      $(HEADERS) $(SRC_DIR)/matrix.cpp
	$(MEX) $(MEXCFLAGS) $(XFLAGS) -I$(INC_DIR) -c $(SRC_DIR)/matrix.cpp

symmatrix.o:   $(HEADERS) $(SRC_DIR)/symmatrix.cpp
	$(MEX) $(MEXCFLAGS) $(XFLAGS) -I$(INC_DIR) -c $(SRC_DIR)/symmatrix.cpp

spmatrix.o:    $(HEADERS) $(SRC_DIR)/spmatrix.cpp
	$(MEX) $(MEXCFLAGS) $(XFLAGS) -I$(INC_DIR) -c $(SRC_DIR)/spmatrix.cpp
endif


# the archieve
OBJS = matkitfunc.o  basics.o vector.o matrix.o symmatrix.o spmatrix.o
$(LIB_DIR)/libsmatkit_mex.a: $(OBJS)
	if [ ! -d $(LIB_DIR) ]; then \
    mkdir $(LIB_DIR); \
fi
	$(AR) $(ARFLAGS) $(LIB_DIR)/libsmatkit_mex.a $(OBJS)

clean:
	rm -f *.o $(LIB_DIR)/libsmatkit_mex.a
