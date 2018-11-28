OBJ = \
rt_tddft_input_put_defaults.o \
rt_tddft_input_field.o \
rt_tddft_input_print.o \
rt_tddft_input_excite.o \
rt_tddft_input_tolerances.o \
rt_tddft_input_visualization.o \
rt_tddft_input_mocap.o \
rt_tddft_input_spatialcap.o \
rt_tddft_input.o


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

