#
# $Id: GNUmakefile,v 1.4 2006-09-18 15:44:09 d3j191 Exp $
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
