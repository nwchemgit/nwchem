#$Id: GNUmakefile,v 1.8 1995-02-02 23:23:10 d3g681 Exp $
# makefile,v 1.4 1994/11/23 18:48:11 og845 Exp

LIBRARY = libutil.a
USES_BLAS = sym_put_in_geom.F sym_print_all.F
HEADERS = sym.fh

OBJ_OPTIMIZE = sym_shell_pair.o 
OBJ =     dosymops.o \
          gensym.o \
          ludcmp.o \
          mprint.o \
          opprint.o \
          sym_put_in_geom.o \
          spgen.o \
          sym_map.o \
          sym_nwc.o \
          sym_apply_op.o \
          sym_get_cart_op.o \
          sym_apply_cart_op.o \
          sym_center_map.o \
          sym_number_ops.o \
          sym_bas_op.o \
          sym_ops_get.o \
          sym_print_all.o \
          sym_symmetrize.o \
          sym_geom_project.o \
          sym_test.o \
          wrcell.o


include ../config/makefile.h
include ../config/makelib.h

