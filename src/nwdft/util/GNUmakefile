# $Id: GNUmakefile,v 1.12 1994-11-24 05:06:03 og845 Exp $

      OBJ_OPTIMIZE = \
	dzero.o		\
	ga_chol.o	\
	ga_inverse.o		\
	gacrea.o	gat_lab_can.o	\
	gather.o	\
	header.o	\
	ijcanon.o	\
	invio.o		\
	pertbl.o	\
	sca_lab_mat.o	scatter.o	\
	second.o	\
	tmpfiles.o	\
	tstamp.o	\
	int_eval_gbsets.o \
	bas_ncontr_cn_max.o		ga_chol_seq.o	\
	ga_inv_seq.o 	\
	int_2e3c_can.o	nbf_to_mbf.o 

     LIBRARY = libnwdft.a

   USES_BLAS = simltr.F int_eval_gbsets.F second.F

include ../../config/makefile.h
include ../../config/makelib.h

ifeq ($(TARGET),SGITFP)
bas_ncontr_cn_max.o:	bas_ncontr_cn_max.F
	f77   $(FOPTIONS) $(INCLUDES) $(DEFINES)  -c bas_ncontr_cn_max.F -o bas_ncontr_cn_max.o
nbf_to_mbf.o:	nbf_to_mbf.F
	f77   $(FOPTIONS) $(INCLUDES) $(DEFINES)  -c nbf_to_mbf.F -o nbf_to_mbf.o
int_eval_gbsets.o:	int_eval_gbsets.F
	f77   $(FOPTIONS) $(INCLUDES) $(DEFINES)  -c int_eval_gbsets.F -o int_eval_gbsets.o
endif
