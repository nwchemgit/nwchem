#
# $Id: GNUmakefile,v 1.3 2006-08-25 20:24:45 d3p307 Exp $
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

ifeq ($(FC),ifc)
FOPTIMIZE += -Ob0
endif

ifeq ($(FC),ifort)
 ifeq ($(_CPU),x86_64)
  FOPTIMIZE += -O0
 endif
endif


include ../config/makelib.h
