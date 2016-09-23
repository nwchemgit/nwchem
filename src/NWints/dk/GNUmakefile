
# $Id$

include ../../config/makefile.h
     OBJ_OPTIMIZE = int_1edk_ga.o \
                    dk_input.o dkdiv_energy.o int_pxvp_ga.o \
                    int_1edkso_ga.o

ifdef XLF10
    OBJ = bas_create_uncont.o
else
    OBJ_OPTIMIZE += bas_create_uncont.o
endif

     LIBRARY = libnwints.a

USES_BLAS = int_1edk_ga.F  int_1edkso_ga.F

include ../../config/makelib.h

