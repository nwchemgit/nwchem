# $Id$

HEADERS = rel_consts.fh \
          rel_nwc.fh \
          zora.fh \
          modelpotential.fh \
          x2c.fh

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
      dmd_input.o \
      zora_input.o \
      x2c_input.o

OBJ_GRAD = \
      rel_oneld.o \
      rel_oneld_cosmo.o \
      rel_2e4cd_sf.o \
      rel_LLSSd.o \
      rel_SSLLd.o \
      rel_SSSSd.o 

OBJ = \
      int_rel_init.o\
      rel_input.o \
      modelpotential_input.o \
      set_modelpotential_params.o \
      get_modelpotential_params.o \
      calc_modelpotential.o

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
      rel_2e4cd_sf.F \
      rel_LLSSd.F \
      rel_SSLLd.F \
      rel_SSSSd.F \
      rel_pot.F  \
      rel_oneld_cosmo.F  \
      rel_pot2.F

LIBRARY = libnwints.a

include ../../config/makefile.h
include ../../config/makelib.h


