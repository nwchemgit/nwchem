# $Id$

HEADERS = libxc.fh

OBJ = nwchem_libxc_read.o \
      nwchem_libxc_util.o
OBJ_OPTIMIZE = nwchem_libxc_compute.o

USES_BLAS = nwchem_libxc_compute.o

LIBRARY = libnwdft.a

LIB_INCLUDES = -I../include 
ifdef USE_LIBXC
LIB_DEFINES += -DXC_MAJOR_VERSION=$(shell grep XC_MAJOR_VERSION ../../libext/libxc/install/include/xc_version.h |cut -d " " -f 3)
LIB_DEFINES += -DLIBXC_NP_SIZE=$(shell ./libxc_findsizenp.sh ../../libext/libxc/install/include)
LIB_INCLUDES += -I../../libext/libxc/install/include
endif
ifdef LIBXC_INCLUDE
LIB_DEFINES += -DXC_MAJOR_VERSION=$(shell grep XC_MAJOR_VERSION $(LIBXC_INCLUDE)/xc_version.h |cut -d " " -f 3)
LIB_DEFINES += -DLIBXC_NP_SIZE=$(shell ./libxc_findsizenp.sh $(LIBXC_INCLUDE))
LIB_INCLUDES += -I$(LIBXC_INCLUDE)
endif
ifdef LIBXC_MODDIR
LIB_INCLUDES += -I$(LIBXC_MODDIR)
endif
LIB_INCLUDES += -I.


include ../../config/makefile.h
include ../../config/makelib.h
