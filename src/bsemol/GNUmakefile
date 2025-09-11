OBJ = bsemol.o bse_input.o bse_init.o bse_analytic.o \
      bse_davidson.o bse_lanczos.o bse_defaults.o \
      bse_finalize.o bse_ri_init.o bse_wmn.o \
      bse_analytic_cis.o bse_analytic_tda.o \
      bse_dipole.o bse_oscstr.o bse_puteris.o bse_get_eia.o \
      bse_davidson_diagonal.o \
      bse_davidson_kvec.o bse_davidson_guess.o \
      bse_davidson_kortho.o bse_davidson_restart.o \
      bse_davidson_newvecs.o bse_davidson_residual.o \
      bse_s2.o

USES_BLAS = bse_buildw.F bse_ri_init.F bse_wmn.F \
            bse_davidson.F bse_davidson_diagonal.F bse_davidson_kvec.F bse_davidson_kortho.F \
            bse_init.F bse_puteris.F bse_s2.F

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
