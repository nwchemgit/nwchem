
# $Id: GNUmakefile,v 1.18 2008-01-25 20:21:07 d3p307 Exp $

include ../../config/makefile.h
     OBJ_OPTIMIZE = bas_create_uncont.o int_1edk_ga.o \
                    dk_input.o dkdiv_energy.o int_pxvp_ga.o \
                    rel_input.o int_1edkso_ga.o

ifdef XLF10
    OBJ = bas_create_uncont.o
else
    OBJ_OPTIMIZE += bas_create_uncont.o
endif

     LIBRARY = libnwints.a

USES_BLAS = int_1edk_ga.F  int_1edkso_ga.F

include ../../config/makelib.h

