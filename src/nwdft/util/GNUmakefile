#$Id: GNUmakefile,v 1.25 1995-05-08 19:46:45 og845 Exp $
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
	bas_ncontr_cn_max.o		ga_chol_seq.o	\
	ga_inv_seq.o 	\
	int_2e3c_can.o	nbf_to_mbf.o dft_mpole.o
endif
else
	OBJ = bas_ncontr_cn_max.o 	nbf_to_mbf.o int_eval_gbsets.o new_eval_gbsets.o
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
	int_2e3c_can.o dft_mpole.o
endif

     LIBRARY = libnwdft.a

   USES_BLAS = int_eval_gbsets.F new_eval_gbsets.F dft_second.F \
               ga_chol_seq.F ga_inv_seq.F dft_mpole.F

include ../../config/makefile.h
include ../../config/makelib.h
