ifndef USE_TBLITE
$(error USE_TBLITE must be set in order to use xTB methods)
endif

ifdef USE_TBLITE
    DEFINES  += -DUSE_TBLITE
    ifdef TBLITE_MESON
        TBLITE_MODS=$(NWCHEM_TOP)/src/libext/tblite/tblite/_build/libtblite.so.0.2.0.p
        MCTC_MODS=$(NWCHEM_TOP)/src/libext/tblite/tblite/_build/subprojects/mctc-lib/libmctc-lib.a.p
    else
        TBLITE_MODS=$(NWCHEM_TOP)/src/libext/include/tblite
        MCTC_MODS=$(sort $(dir $(wildcard $(NWCHEM_TOP)/src/libext/include/mctc-lib/*/)))
        MCHR_MODS=$(sort $(dir $(wildcard $(NWCHEM_TOP)/src/libext/include/multicharge/*/)))
        TOMLF_MODS=$(sort $(dir $(wildcard $(NWCHEM_TOP)/src/libext/include/toml-f/*/)))
        DFTD3_MODS=$(sort $(dir $(wildcard $(NWCHEM_TOP)/src/libext/include/s-dftd3/*/)))
        DFTD4_MODS=$(sort $(dir $(wildcard $(NWCHEM_TOP)/src/libext/include/dftd4/*/)))
    endif
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
