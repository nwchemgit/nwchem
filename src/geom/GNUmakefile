#$Id: GNUmakefile,v 1.8 1995-03-31 01:39:07 d3g681 Exp $

 LIBRARY = libutil.a
#LIBRARY = libgeom.a
     OBJ = geom.o geom_input.o geom_3d.o
 HEADERS = geom.fh geomP.fh

 USES_BLAS = geom.F geom_input.F

include ../config/makefile.h
include ../config/makelib.h

geom_input.o geom.o geom_3d.o:	geomP.fh


