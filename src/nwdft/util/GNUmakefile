# makefile,v 1.15 1994/12/14 02:13:59 og845 Exp
ifeq ($(NWCHEM_TARGET),SGITFP)
	OBJ = bas_ncontr_cn_max.o 	nbf_to_mbf.o \
	int_eval_gbsets.o 
	OBJ_OPTIMIZE = \
	dzero.o		\
	ga_chol.o	\
	ga_inverse.o		\
	gacrea.o	gat_lab_can.o	\
	gather.o	\
	get_col.o	\
	header.o	\
	ijcanon.o	\
	invio.o		\
	pertbl.o	\
	sca_lab_mat.o	scatter.o	\
	second.o	\
	tstamp.o	\
	ga_chol_seq.o	\
	ga_inv_seq.o 	\
	int_2e3c_can.o
else
	OBJ_OPTIMIZE = \
	dzero.o		\
	ga_chol.o	\
	ga_inverse.o		\
	gacrea.o	gat_lab_can.o	\
	gather.o	\
	get_col.o	\
	header.o	\
	ijcanon.o	\
	invio.o		\
	pertbl.o	\
	sca_lab_mat.o	scatter.o	\
	second.o	\
	tstamp.o	\
	int_eval_gbsets.o \
	bas_ncontr_cn_max.o		ga_chol_seq.o	\
	ga_inv_seq.o 	\
	int_2e3c_can.o	nbf_to_mbf.o 
endif
ifeq ($(NWCHEM_TARGET),CRAY-T3D)
	OBJ = ga_chol_seq.o ga_inv_seq.o invio.o
	OBJ_OPTIMIZE = \
	int_eval_gbsets.o \
	dzero.o		\
	ga_chol.o	\
	ga_inverse.o		\
	gacrea.o	gat_lab_can.o	\
	gather.o	\
	get_col.o	\
	header.o	\
	ijcanon.o	\
	pertbl.o	\
	sca_lab_mat.o	scatter.o	\
	second.o	\
	tstamp.o	\
	bas_ncontr_cn_max.o       nbf_to_mbf.o \
	int_2e3c_can.o
endif

     LIBRARY = libnwdft.a

   USES_BLAS = int_eval_gbsets.F second.F ga_chol_seq.F ga_inv_seq.F

include ../../config/makefile.h
include ../../config/makelib.h
