
 LIBRARY = libutil.a
#LIBRARY = libgeom.a
     OBJ = geom.o geom_input.o geom_3d.o
ifeq ($(TARGET),CRAY-T3D)
     OBJ = geom_input.o geom_3d.o
     OBJ_OPTIMIZE = geom.o 
endif
 HEADERS = geom.fh geomP.fh

 USES_BLAS = geom.F geom_input.F

include ../config/makefile.h
include ../config/makelib.h

geom_input.o geom.o geom_3d.o:	geomP.fh


