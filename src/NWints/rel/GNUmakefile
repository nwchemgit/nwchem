# $Id: GNUmakefile,v 1.3 1998-10-23 18:47:08 mg201 Exp $

HEADERS = rel_consts.fh

OBJ_ENERGY = \
      rel_onel.o \
      rel_pot.o  \
      rel_pot2.o \
      rel_LLSS.o \
      rel_SSLL.o \
      rel_SSSS.o \
      rel_init_con.o   

OBJ = 

ifdef NWCHEM_DEBUG
OBJ +=  $(OBJ_ENERGY)

OBJ_OPTIMIZE = 
else

OBJ_OPTIMIZE = $(OBJ_ENERGY)
endif

USES_BLAS = \
      rel_onel.F \
      rel_LLSS.F \
      rel_SSLL.F \
      rel_SSSS.F \
      rel_pot.F  \
      rel_pot2.F

LIBRARY = libnwints.a

include ../../config/makefile.h
include ../../config/makelib.h


