# directories where to "make"
MATKIT_ROOT      = ./MATKIT
FILTLAN_DIRS     = ./OBJ/DOUBLE ./OBJ/SINGLE ./DRIVERS/DOUBLE ./DRIVERS/SINGLE
FILTLAN_MEX_DIRS = ./OBJ_MEX/DOUBLE ./OBJ_MEX/SINGLE ./DRIVERS_MEX


include Make.inc

ifeq ($(MEX_FLAG),NO_MEX)
ALL_DIRS         = $(MATKIT_ROOT) $(FILTLAN_DIRS)
else
ALL_DIRS         = $(MATKIT_ROOT) $(FILTLAN_DIRS) $(FILTLAN_MEX_DIRS)
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
