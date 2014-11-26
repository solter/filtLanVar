# directories where to "make"
MATKIT_ROOT      = .
MATKIT_DIRS      = $(MATKIT_ROOT)/OBJ/DOUBLE  $(MATKIT_ROOT)/OBJ/SINGLE  $(MATKIT_ROOT)/TESTS
MATKIT_MEX_DIRS  = $(MATKIT_ROOT)/OBJ_MEX/DOUBLE  $(MATKIT_ROOT)/OBJ_MEX/SINGLE


include ../Make.inc

ifeq ($(MEX_FLAG),NO_MEX)
ALL_DIRS         = $(MATKIT_DIRS)
else
ALL_DIRS         = $(MATKIT_DIRS) $(MATKIT_MEX_DIRS)
endif


all:
	@for dir in $(ALL_DIRS) ;\
        do\
        echo make all in $$dir;\
        (cd $$dir; make all);\
        done

clean:
	@for dir in $(ALL_DIRS) ;\
        do\
        echo cleaning $$dir ;\
        (cd $$dir;  make clean) ;\
        done
