#$Id: GNUmakefile,v 1.20 1996-01-19 07:32:25 d3g681 Exp $
# makefile,v 1.4 1994/11/23 18:48:11 og845 Exp

LIBRARY = libutil.a
USES_BLAS = sym_put_geom.F sym_pr_all.F sym_g_sym.F \
            sym_mo_adapt.F sym_geom_prj.F sym_bs_irrep.F sym_op_clsfy.F
HEADERS = sym.fh

OBJ_OPTIMIZE = sym_sh_pair.o sym_sym.o \
               sym_mo_ap_op.o sym_bas_op.o
OBJ =     dosymops.o \
          gensym.o \
          ludcmp.o \
          mprint.o \
          opprint.o \
          sym_put_geom.o \
          spgen.o \
          sym_map.o \
          sym_nwc.o \
          sym_apply_op.o \
          sym_get_cart.o \
          sym_grp_name.o \
          sym_ap_cart.o \
          sym_cent_map.o \
          sym_num_ops.o \
          sym_ops_get.o \
          sym_pr_all.o \
          sym_geom_prj.o \
          cross_prod.o \
          sym_op_cname.o \
          sym_pr_ops.o \
          deter3.o \
          sym_op_clsfy.o \
          sym_tr_bs_op.o \
          sym_bs_irrep.o \
          sym_op_type.o \
          sym_char_tab.o \
          sym_pr_ctab.o \
          sym_inv_op.o \
          sym_mo_adapt.o \
          sym_g_sym.o \
          wrcell.o


include ../config/makefile.h
include ../config/makelib.h

ifeq ($(TARGET),CRAY-T3D)
  FNOOPT = -dp -Ccray-t3d -Wf"-o noscalar,jump,noieeedivide"
$(LIBRARY_PATH)(dosymops.o):  dosymops.f
	$(FC) -c $(FNOOPT) $^
$(LIBRARY_PATH)(sym_mo_ap_op.o):  sym_mo_ap_op.F
	$(FCONVERT)
	$(FC) -c $(FNOOPT) sym_mo_ap_op.f
	@/bin/rm -f sym_mo_ap_op.f
endif
             

# DO NOT EDIT BENEATH THIS LINE ... GENERATED AUTOMATICALLY
$(LIBRARY_PATH)(sym_ap_cart.o):	 $(INCDIR)/mafdecls.fh $(INCDIR)/geomP.fh 
$(LIBRARY_PATH)(sym_apply_op.o):	 $(INCDIR)/mafdecls.fh $(INCDIR)/geomP.fh 
$(LIBRARY_PATH)(sym_bas_op.o):	 $(INCDIR)/mafdecls.fh $(INCDIR)/geom.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(sym_bs_irrep.o):	 $(INCDIR)/util.fh ./sym.fh $(INCDIR)/mafdecls.fh $(INCDIR)/geom.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(sym_cent_map.o):	 $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/geomP.fh $(INCDIR)/geom.fh 
$(LIBRARY_PATH)(sym_char_tab.o):	 $(INCDIR)/util.fh 
$(LIBRARY_PATH)(sym_g_sym.o):	 $(INCDIR)/util.fh ./sym.fh $(INCDIR)/mafdecls.fh $(INCDIR)/geom.fh 
$(LIBRARY_PATH)(sym_geom_prj.o):	 $(INCDIR)/util.fh ./sym.fh $(INCDIR)/mafdecls.fh $(INCDIR)/geomP.fh 
$(LIBRARY_PATH)(sym_get_cart.o):	 $(INCDIR)/mafdecls.fh $(INCDIR)/geomP.fh 
$(LIBRARY_PATH)(sym_grp_name.o):	 $(INCDIR)/geomP.fh 
$(LIBRARY_PATH)(sym_inv_op.o):	 ./sym.fh $(INCDIR)/geomP.fh 
$(LIBRARY_PATH)(sym_mo_adapt.o):	 $(INCDIR)/util.fh ./sym.fh $(INCDIR)/msgids.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(sym_mo_ap_op.o):	 ./sym.fh $(INCDIR)/mafdecls.fh $(INCDIR)/geom.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(sym_num_ops.o):	 $(INCDIR)/geomP.fh 
$(LIBRARY_PATH)(sym_nwc.o):	 $(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/geomP.fh $(INCDIR)/geom.fh 
$(LIBRARY_PATH)(sym_op_clsfy.o):	 $(INCDIR)/util.fh ./sym.fh $(INCDIR)/inp.fh 
$(LIBRARY_PATH)(sym_op_cname.o):	 $(INCDIR)/inp.fh 
$(LIBRARY_PATH)(sym_op_mtab.o):	 $(INCDIR)/util.fh ./sym.fh $(INCDIR)/inp.fh 
$(LIBRARY_PATH)(sym_ops_get.o):	 $(INCDIR)/mafdecls.fh $(INCDIR)/geomP.fh 
$(LIBRARY_PATH)(sym_pr_all.o):	 $(INCDIR)/util.fh $(INCDIR)/mafdecls.fh $(INCDIR)/inp.fh $(INCDIR)/geomP.fh 
$(LIBRARY_PATH)(sym_pr_ctab.o):	 ./sym.fh 
$(LIBRARY_PATH)(sym_pr_ops.o):	 ./sym.fh $(INCDIR)/inp.fh 
$(LIBRARY_PATH)(sym_put_geom.o):	 $(INCDIR)/mafdecls.fh $(INCDIR)/geomP.fh $(INCDIR)/geom.fh 
$(LIBRARY_PATH)(sym_sh_pair.o):	 $(INCDIR)/mafdecls.fh $(INCDIR)/itri.fh $(INCDIR)/global.fh $(INCDIR)/geomP.fh $(INCDIR)/geom.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(sym_sym.o):	 $(INCDIR)/tcgmsg.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh $(INCDIR)/bas.fh 
