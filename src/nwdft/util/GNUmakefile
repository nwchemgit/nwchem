#$Id: GNUmakefile,v 1.30 1995-10-25 20:55:44 og845 Exp $
# makefile,v 1.15 1994/12/14 02:13:59 og845 Exp
ifneq ($(NWCHEM_TARGET),SGITFP)
ifeq ($(NWCHEM_TARGET),CRAY-T3D)
	OBJ = ga_chol_seq.o ga_inv_seq.o dft_invio.o
	OBJ_OPTIMIZE = \
	int_eval_gbsets.o 	new_eval_gbsets.o \
	dzero.o		\
	ga_chol.o	\
	ga_inverse.o		\
	gacrea.o	gat_lab_can.o	\
	dft_gather.o	dft_scatter.o \
	get_col.o	\
	dft_header.o	\
	pertbl.o	\
	sca_lab_mat.o	\
	dft_second.o	\
	dft_tstamp.o	\
	bas_ncontr_cn_max.o       nbf_to_mbf.o \
	dft_mpole.o
else
	OBJ_OPTIMIZE = \
	dzero.o		\
	ga_chol.o	\
	ga_inverse.o		\
	gacrea.o	gat_lab_can.o	\
	get_col.o	\
	dft_header.o	\
	dft_invio.o		\
	pertbl.o	\
	dft_gather.o	sca_lab_mat.o	dft_scatter.o	\
	dft_second.o	\
	dft_tstamp.o	\
	int_eval_gbsets.o 	new_eval_gbsets.o \
        dft_eval_basis.o        dft_eval_funs.o   \
        dft_eval_grad.o        dft_eval_hess.o    \
	bas_ncontr_cn_max.o		ga_chol_seq.o	\
	ga_inv_seq.o 	\
	nbf_to_mbf.o dft_mpole.o
endif
else
	OBJ = bas_ncontr_cn_max.o 	nbf_to_mbf.o int_eval_gbsets.o new_eval_gbsets.o         dft_eval_basis.o        dft_eval_funs.o   \
        dft_eval_grad.o        dft_eval_hess.o dft_mpole.o


	OBJ_OPTIMIZE=  \
	dzero.o		\
	ga_chol.o	\
	ga_inverse.o		\
	gacrea.o	gat_lab_can.o	\
	get_col.o	\
	dft_header.o	\
	dft_invio.o		\
	pertbl.o	\
	dft_gather.o	sca_lab_mat.o	dft_scatter.o	\
	dft_second.o	\
	dft_tstamp.o	\
	ga_chol_seq.o	\
	ga_inv_seq.o 	
endif

     LIBRARY = libnwdft.a

   USES_BLAS = int_eval_gbsets.F new_eval_gbsets.F dft_second.F \
               ga_chol_seq.F ga_inv_seq.F dft_mpole.F \
        dft_eval_basis.F        dft_eval_funs.F   \
        dft_eval_grad.F        dft_eval_hess.F 


include ../../config/makefile.h
include ../../config/makelib.h
# DO NOT EDIT BENEATH THIS LINE ... GENERATED AUTOMATICALLY
bas_ncontr_cn_max.o :	 $(INCDIR)/bas.fh 
dft_eval_basis.o :	 $(INCDIR)/mafdecls.fh $(INCDIR)/dimensions.fh $(INCDIR)/bas.fh 
dft_eval_funs.o :	 $(INCDIR)/dimensions.fh $(INCDIR)/bas.fh 
dft_eval_grad.o :	 $(INCDIR)/dimensions.fh $(INCDIR)/bas.fh 
dft_eval_hess.o :	 $(INCDIR)/dimensions.fh $(INCDIR)/bas.fh 
dft_header.o :	 $(INCDIR)/stdio.fh 
dft_invio.o :	 $(INCDIR)/stdio.fh $(INCDIR)/mafdecls.fh 
dft_mpole.o :	 $(INCDIR)/tcgmsg.fh $(INCDIR)/stdio.fh $(INCDIR)/mafdecls.fh $(INCDIR)/inp.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh $(INCDIR)/bas.fh 
dft_second.o :	 $(INCDIR)/timer.fh $(INCDIR)/tcgmsg.fh $(INCDIR)/msgids.fh $(INCDIR)/global.fh 
dft_tstamp.o :	 $(INCDIR)/tcgmsg.fh $(INCDIR)/stdio.fh $(INCDIR)/global.fh 
ga_chol.o :	 $(INCDIR)/tcgmsg.fh $(INCDIR)/msgids.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
ga_chol_seq.o :	 $(INCDIR)/tcgmsg.fh $(INCDIR)/stdio.fh $(INCDIR)/msgids.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
ga_inv_seq.o :	 $(INCDIR)/tcgmsg.fh $(INCDIR)/stdio.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
ga_inverse.o :	 $(INCDIR)/tcgmsg.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
gacrea.o :	 $(INCDIR)/tcgmsg.fh $(INCDIR)/stdio.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
gat_lab_can.o :	 $(INCDIR)/bas.fh 
int_eval_gbsets.o :	 $(INCDIR)/bas.fh 
nbf_to_mbf.o :	 $(INCDIR)/bas.fh 
new_eval_gbsets.o :	 $(INCDIR)/dimensions.fh $(INCDIR)/bas.fh 
sca_lab_mat.o :	 $(INCDIR)/bas.fh 
