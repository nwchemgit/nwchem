# $Id: GNUmakefile,v 1.4 1999-10-26 19:33:41 mg201 Exp $

HEADERS = rel_consts.fh \
          rel_nwc.fh

OBJ_ENERGY = \
      rel_onel.o \
      rel_pvp.o  \
      rel_pot.o  \
      rel_pot2.o \
      rel_LLSS.o \
      rel_SSLL.o \
      rel_SSSS.o \
      rel_init_con.o   

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
      rel_LLSS.F \
      rel_SSLL.F \
      rel_SSSS.F \
      rel_pot.F  \
      rel_pot2.F

LIBRARY = libnwints.a

include ../../config/makefile.h
include ../../config/makelib.h


