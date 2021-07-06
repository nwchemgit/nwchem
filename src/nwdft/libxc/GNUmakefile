# $Id$

HEADERS = libxc.fh

OBJ = nwchem_libxc_read.o \
      nwchem_libxc_util.o
OBJ_OPTIMIZE = nwchem_libxc_compute.o

USES_BLAS = nwchem_libxc_compute.o

LIBRARY = libnwdft.a

LIB_INCLUDES = -I../include -I../../libext/libxc/install/include -I.

include ../../config/makefile.h
include ../../config/makelib.h
