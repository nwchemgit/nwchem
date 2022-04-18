ifndef USE_TBLITE
$(error USE_TBLITE must be set in order to use xTB methods)
endif

OBJ = xtb_energy.o xtb_input.o xtb_energy_gradient.o \
      xtb_rdinput.o xtb_guess.o xtb_spinpol.o xtb_print_param.o \
      xtb_restart_read.o xtb_restart_write.o

LIBRARY = libnwxtb.a

include ../config/makefile.h

LIB_INCLUDES += -I$(TBLITE_MODS)
LIB_INCLUDES += -I$(MCTC_MODS)
LIB_INCLUDES += -I$(MCHR_MODS)
LIB_INCLUDES += -I$(TOMLF_MODS)
LIB_INCLUDES += -I$(DFTD3_MODS)
LIB_INCLUDES += -I$(DFTD4_MODS)
LIB_INCLUDES += -I.

include ../config/makelib.h
