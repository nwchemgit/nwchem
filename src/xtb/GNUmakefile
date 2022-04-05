OBJ = xtb_energy.o xtb_input.o xtb_utils.o xtb_energy_gradient.o \
      xtb_rdinput.o

LIBRARY = libnwxtb.a

include ../config/makefile.h

LIB_INCLUDES += -I$(TBLITE_MODS)
LIB_INCLUDES += -I$(MCTC_MODS)
LIB_INCLUDES += -I.

include ../config/makelib.h
