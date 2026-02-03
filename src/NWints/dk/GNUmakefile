
# $Id$

include ../../config/makefile.h
     OBJ_OPTIMIZE = int_1edk_ga.o \
                    dk_input.o dkdiv_energy.o int_pxvp_ga.o \
                    int_1edkso_ga.o dk_utils.o

ifdef XLF10
    OBJ = bas_create_uncont.o
else
    OBJ_OPTIMIZE += bas_create_uncont.o
endif

     LIBRARY = libnwints.a

USES_BLAS = int_1edk_ga.F  int_1edkso_ga.F

_PDSYEVR := $(shell [ -e ${NWCHEM_TOP}/src/tools/install/lib/libga.a ] &&  nm ${NWCHEM_TOP}/src/tools/install/lib/libga.a|egrep T |egrep  ga_pdsyevr |awk '/ga_pdsyevr/ {print "Y"}')
ifeq ($(_PDSYEVR),Y)
 LIB_DEFINES += -DUSE_PDSYEVR
endif

include ../../config/makelib.h

