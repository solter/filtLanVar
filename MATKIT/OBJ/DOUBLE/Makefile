ROOT_DIR   = ../../..
MATKIT_DIR = ../..

include $(ROOT_DIR)/Make.inc

INC_DIR    = $(MATKIT_DIR)/INCLUDE
SRC_DIR    = $(MATKIT_DIR)/SRC
LIB_DIR    = $(ROOT_DIR)/LIB


all: $(LIB_DIR)/libdmatkit.a


# by default, we use double precision (not single precision) floating point numbers
# XFLAGS being nothing means using the default
XFLAGS =


# the fundamental classes
HEADERS = $(INC_DIR)/matkit.h $(INC_DIR)/matkitfunc.h  $(INC_DIR)/basics.h $(INC_DIR)/vector.h $(INC_DIR)/matrix.h $(INC_DIR)/symmatrix.h $(INC_DIR)/spmatrix.h

matkitfunc.o:  $(INC_DIR)/matkitdef.h $(SRC_DIR)/matkitfunc.cpp
	$(CPP) $(CFLAGS) $(XFLAGS) -I$(INC_DIR) -c $(SRC_DIR)/matkitfunc.cpp

basics.o:      $(INC_DIR)/basics.h $(SRC_DIR)/basics.cpp
	$(CPP) $(CFLAGS) $(XFLAGS) -I$(INC_DIR) -c $(SRC_DIR)/basics.cpp

vector.o:      $(HEADERS) $(SRC_DIR)/vector.cpp
	$(CPP) $(CFLAGS) $(XFLAGS) -I$(INC_DIR) -c $(SRC_DIR)/vector.cpp

matrix.o:      $(HEADERS) $(SRC_DIR)/matrix.cpp
	$(CPP) $(CFLAGS) $(XFLAGS) -I$(INC_DIR) -c $(SRC_DIR)/matrix.cpp

symmatrix.o:   $(HEADERS) $(SRC_DIR)/symmatrix.cpp
	$(CPP) $(CFLAGS) $(XFLAGS) -I$(INC_DIR) -c $(SRC_DIR)/symmatrix.cpp

spmatrix.o:    $(HEADERS) $(SRC_DIR)/spmatrix.cpp
	$(CPP) $(CFLAGS) $(XFLAGS) -I$(INC_DIR) -c $(SRC_DIR)/spmatrix.cpp


# the archieve
OBJS = matkitfunc.o  basics.o vector.o matrix.o symmatrix.o spmatrix.o

$(LIB_DIR)/libdmatkit.a: $(OBJS)
	if [ ! -d $(LIB_DIR) ]; then \
    mkdir $(LIB_DIR); \
fi
	$(AR) $(ARFLAGS) $(LIB_DIR)/libdmatkit.a $(OBJS)


clean:
	rm -f *.o $(LIB_DIR)/libdmatkit.a
