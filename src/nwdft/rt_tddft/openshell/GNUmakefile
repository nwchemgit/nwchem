OBJ = \
zfock_os_build.o \
zfock_os_core.o \
zfock_os_core_scalarzora.o \
zfock_os_coul.o \
zfock_os_coulcdfit.o \
zfock_os_exch.o \
zfock_os_dftxc.o \
rt_tddft_os_init.o \
rt_tddft_os_confirm.o \
rt_tddft_os_prop.o \
rt_tddft_os_tdfock.o \
rt_tddft_os_checkpoint.o \
rt_tddft_os_print.o \
rt_tddft_os_restart.o

OBJ_OPTIMIZE = 

LIBRARY = libnwdft.a

USES_BLAS =

LIB_DEFINES = -DDEBUG_PRINT


##
## KAL: hack to allow symlinked directory
##
LIB_INCLUDES = -I../../include -I../../../ddscf
#LIB_INCLUDES = -I../headers -I$(NWCHEM_TOP)/src/nwdft/include -I$(NWCHEM)/src/ddscf

HEADERS = 


##
## KAL: hack to allow symlinked directory
##
include ../../../config/makefile.h
include ../../../config/makelib.h
#include $(NWCHEM_TOP)/src/config/makefile.h
#include $(NWCHEM_TOP)/src/config/makelib.h
