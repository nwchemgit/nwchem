#$Id: GNUmakefile,v 1.9 1995-04-16 20:21:07 gg502 Exp $

 LIBRARY = libutil.a
#LIBRARY = libgeom.a
     OBJ = geom.o geom_input.o geom_3d.o geom_num_core.o
 HEADERS = geom.fh geomP.fh

 USES_BLAS = geom.F geom_input.F

include ../config/makefile.h
include ../config/makelib.h

geom_input.o geom.o geom_3d.o:	geomP.fh


