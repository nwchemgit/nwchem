#$Id: GNUmakefile,v 1.0 2009-21-06 18:37:00 jhammond Exp $

SUBDIRS = 

OBJ = ncc_input.o ncc_driver.o

OBJ_OPTIMIZE = ncc_util.o ncc_offset_util.o \
               ncc_singles.o ncc_doubles.o \
               ncc_fourindex.o

LIBRARY = libncc.a

USES_BLAS = 

LIB_DEFINES += -DDEBUG_PRINT -DDETAILED_FREE

HEADERS =

include ../config/makefile.h
include ../config/makelib.h

