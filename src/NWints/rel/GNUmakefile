HEADERS   = ecp_nwc.fh

OBJ_ENERGY = \
      rel_onel.o  \
      rel_pot.o   \
      rel_init_con.o   

OBJ = 

ifdef NWCHEM_DEBUG
OBJ +=  $(OBJ_ENERGY)

OBJ_OPTIMIZE = 
else

OBJ_OPTIMIZE = $(OBJ_ENERGY)
endif

HEADERS = ecp_nwc.fh

USES_BLAS = \
      rel_onel.o  \
      rel_pot.o

LIBRARY = libnwints.a

include ../../config/makefile.h
include ../../config/makelib.h


