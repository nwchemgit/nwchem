# $Id: GNUmakefile,v 1.6 1994-11-07 02:09:18 og845 Exp $

      OBJ_OPTIMIZE = \
	dzero.o		\
	ga_chol.o	\
	ga_inverse.o	ga_scac.o	\
	gacrea.o	gat_lab_can.o	\
	gather.o	\
	header.o	\
	ijcanon.o	\
	invio.o		izero.o	\
	pertbl.o	\
	sca_lab_mat.o	sca_lab_can.o		scatter.o	\
	second.o	\
	tmpfiles.o	\
	tstamp.o	\
	int_eval_gbsets.o \
	bas_ncontr_cn_max.o		ga_chol_seq.o	\
	ga_inv_seq.o \
	int_2e3c_can.o	nbf_to_mbf.o

     LIBRARY = libnwdft.a

#   USES_BLAS = 

include ../../config/makefile.h
include ../../config/makelib.h

ifeq ($(TARGET),SGITFP)
bas_ncoef_cn_max.o:	bas_ncoef_cn_max.F
	f77   $(FOPTIONS) $(INCLUDES) $(DEFINES)  -c bas_ncoef_cn_max.F -o bas_ncoef_cn_max.o
endif
