#
# $Id: GNUmakefile,v 1.2 2003-08-11 19:43:18 d3j191 Exp $
#
#
OBJ_OPTIMIZE=qhop_setup.o qhop_asite.o qhop_dsite.o qhop_prob.o qhop_init.o qhop_pindex.o qhop_final.o

LIBRARY=libqhop.a

HEADER=

USES_BLAS=

include ../config/makefile.h

ifeq ($(TARGET),LINUX)
DEFINES += -DNEED_LOC
endif

include ../config/makelib.h
