# $Id: GNUmakefile,v 1.58 2000-01-15 02:05:19 edo Exp $
ifeq ($(TARGET),SGITFP)
	OBJ = xc_getv.o xc_quadv0.o xc_quadv0_a.o xc_quadv1.o xc_fitv.o xc_get2oi.o qpts.o gat_mat.o	xc_beckew.o setACmat.o xc_neighb.o xc_ssw.o xc_oep.o xc_sic_pert.o xc_sic.o localize_sic.o   xc_exact_pot.o
	OBJ_OPTIMIZE = xc_tabcd.o  xc_rhogen.o\
	xc_becke88.o	xc_perdew86.o xc_perdew91.o xc_pw91lda.o\
	xc_dirac.o	xc_vwn_5.o	xc_lyp88.o			\
	eumac.o	\
	scat_mat.o	 \
	qpts.o		   \
	a_dist.o	qdist.o		 xc_delleyw.o xc_lebedev.o \
        xc_vwn_1_rpa.o xc_vwn_1.o xc_vwn_2.o xc_vwn_3.o xc_vwn_4.o \
	xc_eval_fnl.o xc_pbe96.o xc_hcth.o xc_chkgrad.o xc_xpbe96.o xc_xpw91.o xc_tarad.o xc_murakn.o xc_signf.o xc_interp.o xc_gill96.o xc_oep.o xc_sic_pert.o xc_sic.o localize_sic.o 
else
#      OBJ = xc_getv.o  xc_quadv1.o  qpts.o  xc_quadv0.o xc_quadv0_a.o xc_perdew91.o 
#
#
#      OBJ_OPTIMIZE = xc_rhogen.o		xc_pw91lda.o   	\
#	a_dist.o	xc_beckew.o 	 xc_neighb.o\
#	xc_perdew86.o	xc_dirac.o	\
#	xc_becke88.o	xc_lyp88.o	xc_vwn_5.o \
#	eumac.o	\
#	xc_fitv.o			   xc_get2oi.o	\
#	xc_lebedev.o		   \
#	qdist.o	\
#	setACmat.o		xc_tabcd.o	\
#	gat_mat.o		   scat_mat.o xc_delleyw.o \
#        xc_vwn_1_rpa.o xc_vwn_1.o xc_vwn_2.o xc_vwn_3.o xc_vwn_4.o
#
	OBJ=xc_getv.o xc_exact_pot.o
	OBJ_OPTIMIZE = \
	a_dist.o \
	eumac.o \
	gat_mat.o \
	qdist.o \
	qpts.o \
	scat_mat.o \
	setACmat.o \
	xc_becke88.o \
	xc_beckew.o \
	xc_delleyw.o \
	xc_dirac.o \
	xc_fitv.o \
	xc_get2oi.o \
	xc_lebedev.o \
	xc_lyp88.o \
	xc_neighb.o \
	xc_perdew86.o \
	xc_perdew91.o \
	xc_pw91lda.o \
	xc_quadv0.o \
	xc_quadv0_a.o \
	xc_quadv1.o \
	xc_rhogen.o \
	xc_tabcd.o \
	xc_vwn_1.o \
	xc_vwn_1_rpa.o \
	xc_vwn_2.o \
	xc_vwn_3.o \
	xc_vwn_4.o \
	xc_vwn_5.o \
	xc_eval_fnl.o xc_pbe96.o   xc_chkgrad.o\
	xc_hcth.o  xc_xpbe96.o xc_ssw.o  xc_xpw91.o xc_tarad.o xc_murakn.o xc_signf.o xc_interp.o xc_gill96.o xc_oep.o xc_sic_pert.o xc_sic.o localize_sic.o 

endif

     LIBRARY = libnwdft.a

 LIB_DEFINES = -DDFT_TIMER
# LIB_DEFINES = -DQUADV0_REVERT

LIB_INCLUDES = -I../include


   USES_BLAS = xc_fitv.F xc_quadv0_a.F xc_quadv1.F xc_dirac.F xc_tabcd.F xc_rhogen.F xc_beckew.F xc_ssw.F xc_oep.F xc_sic_pert.F xc_sic.F localize_sic.F \
               ../include/xc_rhogen.fh ../include/xc_tabcd.fh

include ../../config/makefile.h
include ../../config/makelib.h


