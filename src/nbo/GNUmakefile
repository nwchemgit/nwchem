#
# $Id: GNUmakefile,v 1.3 2004-04-29 23:01:54 edo Exp $
#


   OBJ = nwnbo.o nbo_holler.o

   LIBRARY = libnbo.a
   LIB_DEFINES = 
   LIB_INCLUDES = -I../include

USES_BLAS = nwnbo.F 

include ../config/makefile.h
include ../config/makelib.h

ifeq ($(_IFCV8),Y)
 FOPTIMIZE += -O0 -g
 FDEBUG += -O0 -g
endif


