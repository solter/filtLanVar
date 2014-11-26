#-----------------------------------------------------------------------------
# platform
#-----------------------------------------------------------------------------

PLAT = LINUX


#-----------------------------------------------------------------------------
# ARPACK (for Lanczos w/ implicit restarting and w/ full reorthogonalization)
#-----------------------------------------------------------------------------

ARPACK_LIB = -larpack
# If you compiled ARPACK with an archive file by yourself, specify it such as
# follows
#ARPACK_DIR = ../../ARPACK
#ARPACK_LIB = $(ARPACK_DIR)/libarpack_$(PLAT).a



#-----------------------------------------------------------------------------
# UMFPACK and supporting libraries (used by the shift-and-invert Lanczos)
#-----------------------------------------------------------------------------

DSDIR  = ../../DIRECT_SPARSE

UFCONFIG_INC = $(DSDIR)/UFconfig
UMF_INC = $(DSDIR)/UMFPACK/Include
AMD_INC = $(DSDIR)/AMD/Include
UMFPACK_INCX = -I$(UFCONFIG_INC) -I$(UMF_INC) -I$(AMD_INC)

UMF_LIB = $(DSDIR)/UMFPACK/Lib/libumfpack.a
AMD_LIB = $(DSDIR)/AMD/Lib/libamd.a
CHOLMOD_LIB = $(DSDIR)/CHOLMOD/Lib/libcholmod.a
CAMD_LIB = $(DSDIR)/CAMD/Lib/libcamd.a
CCOLAMD_LIB = $(DSDIR)/CCOLAMD/Lib/libccolamd.a
COLAMD_LIB = $(DSDIR)/COLAMD/Lib/libcolamd.a
METIS_LIB = $(DSDIR)/metis-4.0/libmetis.a
UMFPACK_LIBX = $(UMF_LIB) $(AMD_LIB) $(CHOLMOD_LIB) $(CAMD_LIB) $(CCOLAMD_LIB) $(COLAMD_LIB) $(METIS_LIB)
