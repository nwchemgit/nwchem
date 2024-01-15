# $Id: GNUmakefile 25498 2014-04-18 02:44:47Z d3y133 $
  LIBRARY = libmm.a

  OBJ =  mm_input.o mm_data.o mm_coords.o mm_vdw.o mm_bond.o mm_bond_coords.o \
         mm_link.o mm_geom.o mm_vdw_coords.o mm_bq.o mm_neb.o mm_cons.o

  USES_BLAS = mm_input.F mm_bond.F mm_bond_coords.F mm_data.F mm_geom.F mm_link.F mm_vdw.F mm_vdw_coords.F

HEADERS = mm.fh mm_params.fh

include ../config/makefile.h
include ../config/makelib.h

