# $Id: GNUmakefile,v 1.6 2000-01-25 19:44:22 mg201 Exp $

HEADERS = rel_consts.fh \
          rel_nwc.fh

OBJ_ENERGY = \
      rel_onel.o \
      rel_mpole.o \
      rel_pvp.o  \
      rel_pot.o  \
      rel_pot2.o \
      rel_2e4c_sf.o \
      rel_LLSS.o \
      rel_SSLL.o \
      rel_SSSS.o \
      rel_init_con.o \
      dmd_input.o

OBJ_GRAD = \
      rel_oneld.o 

OBJ = \
      int_rel_init.o

ifdef NWCHEM_DEBUG
OBJ +=  $(OBJ_ENERGY) $(OBJ_GRAD)

OBJ_OPTIMIZE = 
else

OBJ_OPTIMIZE = $(OBJ_ENERGY) $(OBJ_GRAD)
endif

USES_BLAS = \
      rel_onel.F \
      rel_oneld.F \
      rel_mpole.F \
      rel_2e4c_sf.F \
      rel_LLSS.F \
      rel_SSLL.F \
      rel_SSSS.F \
      rel_pot.F  \
      rel_pot2.F

LIBRARY = libnwints.a

include ../../config/makefile.h
include ../../config/makelib.h


