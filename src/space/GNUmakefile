#
# $Id: GNUmakefile,v 1.5 2007-04-18 18:09:09 d3p307 Exp $
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

ifeq ($(TARGET),LINUX64)
ifeq ($(_CPU),x86_64)
FOPTIMIZE+=-O1
endif
endif

include ../config/makelib.h
