OBJ = \
zfock_cs_build.o \
zfock_cs_core.o \
zfock_cs_core_scalarzora.o \
zfock_cs_coulcdfit.o \
zfock_cs_coul_exchre.o \
zfock_cs_coul.o \
zfock_cs_dftxc.o \
zfock_cs_exch.o \
zfock_cs_exchim.o \
rt_tddft_cs_init.o \
rt_tddft_cs_confirm.o \
rt_tddft_cs_checkpoint.o \
rt_tddft_cs_restart.o \
rt_tddft_cs_prop.o \
rt_tddft_cs_tdfock.o



OBJ_OPTIMIZE = 

LIBRARY = libnwdft.a

USES_BLAS =

LIB_DEFINES = -DDEBUG_PRINT


##
## KAL: hack to allow symlinked directory
##
LIB_INCLUDES = -I../../include -I../../../ddscf
#LIB_INCLUDES = -I../headers -I$(NWCHEM_TOP)/src/nwdft/include -I$(NWCHEM)/src/ddscf

#xxx for cgridfile
#LIB_INCLUDES = -I../headers -I$(NWCHEM_TOP)/src/nwdft/include -I$(NWCHEM_TOP)/src/ddscf -I$(NWCHEM_TOP)/src/nwdft/grid

HEADERS = 


##
## KAL: hack to allow symlinked directory
##
include ../../../config/makefile.h
include ../../../config/makelib.h
#include $(NWCHEM_TOP)/src/config/makefile.h
#include $(NWCHEM_TOP)/src/config/makelib.h
