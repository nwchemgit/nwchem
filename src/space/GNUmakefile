#
# $Id: GNUmakefile,v 1.2 1999-07-27 21:09:07 d3e129 Exp $
#
#
OBJ_OPTIMIZE=sp_init.o sp_main.o

LIBRARY=libspace.a

HEADER=

USES_BLAS=

include ../config/makefile.h

ifeq ($(TARGET),LINUX)
DEFINES += -DNEED_LOC
endif

include ../config/makelib.h
