#
# $Id: GNUmakefile,v 1.1 2003-07-17 19:36:51 d3j191 Exp $
#
#
OBJ_OPTIMIZE=qhop_setup.o qhop_asite.o qhop_dsite.o qhop_prob.o qhop_init.o qhop_pindex.o

LIBRARY=libqhop.a

HEADER=

USES_BLAS=

include ../config/makefile.h

ifeq ($(TARGET),LINUX)
DEFINES += -DNEED_LOC
endif

include ../config/makelib.h
