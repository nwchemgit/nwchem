# $Id: GNUmakefile,v 1.30 1996-01-23 00:29:27 d3e129 Exp $
ifeq ($(TARGET),SGITFP)
	OBJ = xc_getv.o xc_quadv0.o xc_quadv0_a.o xc_quadv1.o xc_fitv.o xc_get2oi.o qpts.o gat_mat.o	xc_beckew.o setACmat.o xc_neighb.o
	OBJ_OPTIMIZE = xc_tabcd.o  xc_rhogen.o\
	xc_becke88.o	xc_perdew86.o xc_perdew91.o \
	xc_dirac.o	xc_vwn.o	xc_lyp88.o			\
	eumac.o	\
	scat_mat.o	 \
	qpts.o		gausquad.o			   \
	a_dist.o	qdist.o		 xc_delleyw.o xc_lebedev.o
else
      OBJ = xc_getv.o  xc_quadv1.o  qpts.o  xc_quadv0.o xc_quadv0_a.o xc_perdew91.o 	


      OBJ_OPTIMIZE = xc_rhogen.o	   	\
	a_dist.o	  xc_beckew.o 	 xc_neighb.o\
	xc_perdew86.o	 xc_dirac.o	\
	xc_becke88.o 	xc_lyp88.o	xc_vwn.o \
	eumac.o	\
	xc_fitv.o			   xc_get2oi.o	\
	gausquad.o	xc_lebedev.o		   \
	qdist.o	\
	setACmat.o		xc_tabcd.o	\
	gat_mat.o		   scat_mat.o xc_delleyw.o
endif

     LIBRARY = libnwdft.a
 LIB_DEFINES = -DDFT_TIMER
LIB_INCLUDES = -I../include


   USES_BLAS = xc_fitv.F xc_quadv0_a.F xc_quadv1.F xc_dirac.F xc_tabcd.F xc_rhogen.F

include ../../config/makefile.h
include ../../config/makelib.h


