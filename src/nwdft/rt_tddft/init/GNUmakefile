OBJ = \
rt_tddft_init_coulcdfit.o \
rt_tddft_init_zora.o \
rt_tddft_init_dip_quad_ints.o \
rt_tddft_init_fock2e.o \
rt_tddft_init_geoms.o \
rt_tddft_init_rtdb.o \
rt_tddft_init_check_params.o \
rt_tddft_init_overlap_canorg.o \
rt_tddft_init_common.o \
rt_tddft_clean.o \
task_rt_tddft.o 

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

