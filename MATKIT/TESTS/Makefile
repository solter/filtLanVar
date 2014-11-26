ROOT_DIR = ../..

include $(ROOT_DIR)/Make.inc

INC_DIR = ../INCLUDE
LIB_DIR = $(ROOT_DIR)/LIB
HEADERS = $(INC_DIR)/matkitdef.h $(INC_DIR)/matkitfunc.h  $(INC_DIR)/basics.h $(INC_DIR)/vector.h $(INC_DIR)/matrix.h $(INC_DIR)/symmatrix.h $(INC_DIR)/spmatrix.h

DFLAGS =
SFLAGS = -DUSE_SINGLE


all: test_matkit test_matkit_efficiency test_matkit_single test_matkit_efficiency_single


# double precision test drivers
test_matkit: test_matkit.o $(LIB_DIR)/libdmatkit.a
	$(CPP) -o test_matkit test_matkit.o $(LIB_DIR)/libdmatkit.a $(LOADFLAGS)

test_matkit.o: test_matkit.cpp $(HEADERS)
	$(CPP) $(CFLAGS) $(DFLAGS) -I$(INC_DIR) -c -o test_matkit.o test_matkit.cpp

test_matkit_efficiency: test_matkit_efficiency.o $(LIB_DIR)/libdmatkit.a
	$(CPP) -o test_matkit_efficiency test_matkit_efficiency.o $(LIB_DIR)/libdmatkit.a $(LOADFLAGS)

test_matkit_efficiency.o: test_matkit_efficiency.cpp $(HEADERS)
	$(CPP) $(CFLAGS) $(DFLAGS) -I$(INC_DIR) -c -o test_matkit_efficiency.o test_matkit_efficiency.cpp


# single precision test drivers
test_matkit_single: test_matkit_single.o $(LIB_DIR)/libsmatkit.a
	$(CPP) -o test_matkit_single test_matkit_single.o $(LIB_DIR)/libsmatkit.a $(LOADFLAGS)

test_matkit_single.o: test_matkit.cpp $(HEADERS)
	$(CPP) $(CFLAGS) $(SFLAGS) -I$(INC_DIR) -c -o test_matkit_single.o test_matkit.cpp

test_matkit_efficiency_single: test_matkit_efficiency_single.o $(LIB_DIR)/libsmatkit.a
	$(CPP) -o test_matkit_efficiency_single test_matkit_efficiency_single.o $(LIB_DIR)/libsmatkit.a $(LOADFLAGS)

test_matkit_efficiency_single.o: test_matkit_efficiency.cpp $(HEADERS)
	$(CPP) $(CFLAGS) $(SFLAGS) -I$(INC_DIR) -c -o test_matkit_efficiency_single.o test_matkit_efficiency.cpp


clean:
	rm -f *.o *~ test_matkit test_matkit_efficiency test_matkit_single test_matkit_efficiency_single
