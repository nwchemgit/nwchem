#$Id: GNUmakefile,v 1.12 1995-10-17 08:56:20 d3g681 Exp $
# makefile,v 1.4 1994/11/23 18:48:11 og845 Exp

LIBRARY = libutil.a
USES_BLAS = sym_put_in_geom.F sym_print_all.F sym_grad_symmetrize.F \
            sym_movecs_adapt.F sym_geom_project.F sym_bas_irreps.F
HEADERS = sym.fh

OBJ_OPTIMIZE = sym_shell_pair.o sym_symmetrize.o \
               sym_movecs_apply_op.o
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
          sym_group_name.o \
          sym_apply_cart_op.o \
          sym_center_map.o \
          sym_number_ops.o \
          sym_bas_op.o \
          sym_ops_get.o \
          sym_print_all.o \
          sym_geom_project.o \
          cross_product.o \
          sym_op_class_name.o \
          sym_print_ops.o \
          deter3.o \
          sym_op_classify.o \
          sym_trace_bas_op.o \
          sym_bas_irreps.o \
          sym_op_type.o \
          sym_char_table.o \
          sym_print_char_table.o \
          sym_inv_op.o \
          sym_movecs_adapt.o \
          sym_grad_symmetrize.o \
          wrcell.o


include ../config/makefile.h
include ../config/makelib.h

# DO NOT EDIT BENEATH THIS LINE ... GENERATED AUTOMATICALLY
sym_apply_cart_op.o :	 $(INCDIR)/mafdecls.fh $(INCDIR)/geomP.fh 
sym_apply_op.o :	 $(INCDIR)/mafdecls.fh $(INCDIR)/geomP.fh 
sym_bas_irreps.o :	 $(INCDIR)/util.fh sym.fh $(INCDIR)/mafdecls.fh $(INCDIR)/geom.fh $(INCDIR)/bas.fh 
sym_bas_op.o :	 $(INCDIR)/mafdecls.fh $(INCDIR)/geom.fh $(INCDIR)/bas.fh 
sym_center_map.o :	 $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/geomP.fh $(INCDIR)/geom.fh 
sym_geom_project.o :	 $(INCDIR)/util.fh sym.fh $(INCDIR)/mafdecls.fh $(INCDIR)/geomP.fh 
sym_get_cart_op.o :	 $(INCDIR)/mafdecls.fh $(INCDIR)/geomP.fh 
sym_grad_symmetrize.o :	 $(INCDIR)/util.fh sym.fh $(INCDIR)/mafdecls.fh $(INCDIR)/geom.fh 
sym_group_name.o :	 $(INCDIR)/geomP.fh 
sym_inv_op.o :	 sym.fh 
sym_movecs_adapt.o :	 $(INCDIR)/util.fh sym.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh $(INCDIR)/bas.fh 
sym_movecs_apply_op.o :	 sym.fh $(INCDIR)/mafdecls.fh $(INCDIR)/geom.fh $(INCDIR)/bas.fh 
sym_number_ops.o :	 $(INCDIR)/geomP.fh 
sym_nwc.o :	 $(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/geomP.fh $(INCDIR)/geom.fh 
sym_op_class_name.o :	 $(INCDIR)/inp.fh 
sym_op_classify.o :	 sym.fh $(INCDIR)/inp.fh 
sym_ops_get.o :	 $(INCDIR)/mafdecls.fh $(INCDIR)/geomP.fh 
sym_print_all.o :	 $(INCDIR)/util.fh $(INCDIR)/mafdecls.fh $(INCDIR)/inp.fh $(INCDIR)/geomP.fh 
sym_print_char_table.o :	 sym.fh 
sym_print_ops.o :	 sym.fh $(INCDIR)/inp.fh 
sym_put_in_geom.o :	 $(INCDIR)/mafdecls.fh $(INCDIR)/geomP.fh $(INCDIR)/geom.fh 
sym_shell_pair.o :	 $(INCDIR)/mafdecls.fh $(INCDIR)/itri.fh $(INCDIR)/global.fh $(INCDIR)/geomP.fh $(INCDIR)/geom.fh $(INCDIR)/bas.fh 
sym_symmetrize.o :	 $(INCDIR)/tcgmsg.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh $(INCDIR)/bas.fh 
