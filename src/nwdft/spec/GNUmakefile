# $Id$
#
	OBJ_OPTIMIZE = \
	dft_vspec.o\
	dft_vspec_so.o\
	dft_xps.o

     LIBRARY = libnwdft.a
_PDSYEVR = $(shell cat ${NWCHEM_TOP}/src/tools/build/config.h | awk ' /HAVE_PDSYEVR\ 1/ {print "Y"}')
ifeq ($(_PDSYEVR),Y) 
 LIB_DEFINES += -DUSE_PDSYVER
endif

LIB_INCLUDES = -I../include -I../grid

USES_BLAS = 

HEADERS =

include ../../config/makefile.h
include ../../config/makelib.h
