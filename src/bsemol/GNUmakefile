OBJ = bsemol.o bse_input.o bse_init.o bse_analytic.o \
      bse_davidson.o bse_lanczos.o bse_defaults.o \
      bse_finalize.o bse_ri_init.o bse_wmn.o \
      bse_analytic_cis.o bse_analytic_tda.o \
      bse_dipole.o bse_oscstr.o bse_puteris.o bse_get_eia.o

USES_BLAS = bse_buildw.F bse_ri_init.F bse_wmn.F

OBJ_OPTIMIZE = bse_buildw.o

LIBRARY = libbsemol.a

HEADERS = ../nwdft/include/cdft.fh ../nwdft/grid/cgridfile.fh ../nwdft/include/dft_scaleig.fh

LIB_DEFINES =
LIB_INCLUDES = -I./

include ../config/makefile.h


ifeq ($(_USE_SCALAPACK),Y)
_USE_ELPA = $(shell ${GA_PATH}/bin/ga-config  --use_elpa| awk ' /1/ {print "Y"}')

ifeq ($(_USE_ELPA),Y)
 LIB_DEFINES += -DUSE_ELPA
endif
#assuming we have GAPZHEEVD
GOTGAPZHEEVD= Y
ifeq ($(GOTGAPZHEEVD),Y)
 LIB_DEFINES += -DUSE_GAPZHEEVD
endif
_PDSYEVR = $(shell [ -e ${NWCHEM_TOP}/src/tools/install/lib/libga.a ] &&  nm ${NWCHEM_TOP}/src/tools/install/lib/libga.a|egrep T |egrep  ga_pdsyevr |awk '/ga_pdsyevr/ {print "Y"}')
ifeq ($(_PDSYEVR),Y) 
 LIB_DEFINES += -DUSE_PDSYEVR
endif
endif
include ../config/makelib.h
