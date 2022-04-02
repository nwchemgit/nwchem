OBJ = xtb_energy.o xtb_input.o

LIBRARY = libnwxtb.a

include ../config/makefile.h

ifdef TBLITE_MODS
LIB_INCLUDES += -I$(TBLITE_MODS)
endif
ifdef MCTC_MODS
LIB_INCLUDES += -I$(MCTC_MODS)
endif

LIB_INCLUDES += -I.

include ../config/makelib.h
