ROOT_DIR      = ../..
LIB_DIR       = $(ROOT_DIR)/LIB

include $(ROOT_DIR)/Make.inc

MATKIT_DIR    = $(ROOT_DIR)/MATKIT
MATKIT_INC    = $(MATKIT_DIR)/INCLUDE
MATKIT_HEADS  = $(MATKIT_INC)/matkitdef.h $(MATKIT_INC)/matkitfunc.h  $(MATKIT_INC)/basics.h $(MATKIT_INC)/vector.h $(MATKIT_INC)/matrix.h $(MATKIT_INC)/symmatrix.h $(MATKIT_INC)/spmatrix.h
MATKIT_LIB    = $(LIB_DIR)/libsmatkit.a

FILTLAN_DIR   = $(ROOT_DIR)
FILTLAN_INC   = $(FILTLAN_DIR)/INCLUDE
LANCZOS_HEADS = $(FILTLAN_INC)/laneig.h $(FILTLAN_INC)/symtrieig.h
FILTLAN_HEADS = $(FILTLAN_INC)/filtlan.h $(FILTLAN_INC)/polyfilt.h $(FILTLAN_INC)/laneig.h $(FILTLAN_INC)/symtrieig.h
FILTLAN_LIB   = $(LIB_DIR)/libsfiltlan.a

DRIVERS_SRC   = $(FILTLAN_DIR)/DRIVERS


DRIVERS = laneig_driver filtlan_driver polyfilt_driver numeig_driver
#DRIVERS = laneig_driver filtlan_driver numeig_driver

all: $(DRIVERS)

XFLAGS = -DUSE_SINGLE

laneig_driver:      laneig_driver.o $(FILTLAN_LIB) $(MATKIT_LIB)
	$(CPP) -o laneig_driver laneig_driver.o $(FILTLAN_LIB) $(MATKIT_LIB) $(LOADFLAGS)

laneig_driver.o:    $(DRIVERS_SRC)/laneig_driver.cpp $(MATKIT_HEADS) $(FILTLAN_INC)/laneig.h $(FILTLAN_INC)/symtrieig.h
	$(CPP) $(CFLAGS) $(XFLAGS) -I$(MATKIT_INC) -I$(FILTLAN_INC) -c $(DRIVERS_SRC)/laneig_driver.cpp


polyfilt_driver:    polyfilt_driver.o $(FILTLAN_LIB) $(MATKIT_LIB)
	$(CPP) -o polyfilt_driver polyfilt_driver.o $(FILTLAN_LIB) $(MATKIT_LIB) $(LOADFLAGS)

polyfilt_driver.o:  $(DRIVERS_SRC)/polyfilt_driver.cpp $(MATKIT_HEADS) $(FILTLAN_INC)/polyfilt.h
	$(CPP) $(CFLAGS) $(XFLAGS) -I$(MATKIT_INC) -I$(FILTLAN_INC) -c $(DRIVERS_SRC)/polyfilt_driver.cpp


filtlan_driver:     filtlan_driver.o $(FILTLAN_LIB) $(MATKIT_LIB)
	$(CPP) -o filtlan_driver filtlan_driver.o $(FILTLAN_LIB) $(MATKIT_LIB) $(LOADFLAGS)

filtlan_driver.o:   $(DRIVERS_SRC)/filtlan_driver.cpp $(MATKIT_HEADS) $(FILTLAN_HEADS)
	$(CPP) $(CFLAGS) $(XFLAGS) -I$(MATKIT_INC) -I$(FILTLAN_INC) -c $(DRIVERS_SRC)/filtlan_driver.cpp


numeig_driver:      numeig_driver.o $(FILTLAN_LIB) $(MATKIT_LIB)
	$(CPP) -o numeig_driver numeig_driver.o $(FILTLAN_LIB) $(MATKIT_LIB) $(LOADFLAGS)

numeig_driver.o:    $(DRIVERS_SRC)/numeig_driver.cpp $(MATKIT_HEADS) $(FILTLAN_HEADS)
	$(CPP) $(CFLAGS) $(XFLAGS) -I$(MATKIT_INC) -I$(FILTLAN_INC) -c $(DRIVERS_SRC)/numeig_driver.cpp


clean:
	rm -f *.o $(DRIVERS)
