# $Id$

HEADERS = libxc.fh

OBJ = nwchem_libxc_read.o \
      nwchem_libxc_util.o
OBJ_OPTIMIZE = nwchem_libxc_compute.o

USES_BLAS = nwchem_libxc_compute.o

LIBRARY = libnwdft.a

LIB_INCLUDES = -I../include 
ifdef USE_LIBXC
LIB_INCLUDES += -I../../libext/libxc/install/include
endif
ifdef LIBXC_DIR
LIB_INCLUDES += -I$(LIBXC_DIR)/include
endif
LIB_INCLUDES += -I.

include ../../config/makefile.h
include ../../config/makelib.h
