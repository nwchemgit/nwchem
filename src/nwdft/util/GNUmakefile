#$Id: GNUmakefile,v 1.28 1995-08-25 19:23:05 d3e129 Exp $
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
	int_2e3c_can.o dft_mpole.o
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
	int_2e3c_can.o	nbf_to_mbf.o dft_mpole.o
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
	ga_inv_seq.o 	\
	int_2e3c_can.o 
endif

     LIBRARY = libnwdft.a

   USES_BLAS = int_eval_gbsets.F new_eval_gbsets.F dft_second.F \
               ga_chol_seq.F ga_inv_seq.F dft_mpole.F \
        dft_eval_basis.F        dft_eval_funs.F   \
        dft_eval_grad.F        dft_eval_hess.F 


include ../../config/makefile.h
include ../../config/makelib.h
