OBJ = bsemol.o bse_input.o

USES_BLAS =
OBJ_OPTIMIZE =

LIBRARY = libbsemol.a

HEADERS = ../nwdft/include/cdft.fh ../nwdft/grid/cgridfile.fh 

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
