# $Id: GNUmakefile,v 1.52 1999-05-10 22:41:52 edo Exp $
ifeq ($(TARGET),SGITFP)
	OBJ = xc_getv.o xc_quadv0.o xc_quadv0_a.o xc_quadv1.o xc_fitv.o xc_get2oi.o qpts.o gat_mat.o	xc_beckew.o setACmat.o xc_neighb.o xc_ssw.o
	OBJ_OPTIMIZE = xc_tabcd.o  xc_rhogen.o\
	xc_becke88.o	xc_perdew86.o xc_perdew91.o xc_pw91lda.o\
	xc_dirac.o	xc_vwn_5.o	xc_lyp88.o			\
	eumac.o	\
	scat_mat.o	 \
	qpts.o		gausquad.o			   \
	a_dist.o	qdist.o		 xc_delleyw.o xc_lebedev.o \
        xc_vwn_1_rpa.o xc_vwn_1.o xc_vwn_2.o xc_vwn_3.o xc_vwn_4.o \
	xc_eval_fnl.o xc_pbe96.o xc_hcth.o xc_chkgrad.o xc_xpbe96.o xc_xpw91.o xc_tarad.o xc_murakn.o xc_signf.o   xc_interp.o xc_octant.o
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
#	gausquad.o	xc_lebedev.o		   \
#	qdist.o	\
#	setACmat.o		xc_tabcd.o	\
#	gat_mat.o		   scat_mat.o xc_delleyw.o \
#        xc_vwn_1_rpa.o xc_vwn_1.o xc_vwn_2.o xc_vwn_3.o xc_vwn_4.o
#

	OBJ_OPTIMIZE = \
	a_dist.o \
	eumac.o \
	gat_mat.o \
	gausquad.o \
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
	xc_getv.o \
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
	xc_hcth.o  xc_xpbe96.o xc_ssw.o  xc_xpw91.o xc_tarad.o xc_murakn.o xc_signf.o   xc_interp.o 
	OBJ=xc_octant.o

endif

     LIBRARY = libnwdft.a

# LIB_DEFINES = -DDFT_TIMER
# LIB_DEFINES = -DQUADV0_REVERT

LIB_INCLUDES = -I../include


   USES_BLAS = xc_fitv.F xc_quadv0_a.F xc_quadv1.F xc_dirac.F xc_tabcd.F xc_rhogen.F xc_beckew.F xc_ssw.F\
               ../include/xc_rhogen.fh ../include/xc_tabcd.fh

include ../../config/makefile.h
include ../../config/makelib.h


