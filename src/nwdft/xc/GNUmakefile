# $Id: GNUmakefile,v 1.20 1995-05-19 00:22:25 og845 Exp $
ifeq ($(NWCHEM_TARGET),SGITFP)
	OBJ = xc_getv.o xc_quadv0.o xc_quadv0_a.o xc_quadv1.o xc_fitv.o xc_get2oi.o qpts.o gat_mat.o	xc_beckew.o
	OBJ_OPTIMIZE = \
	xc_becke88.o	xc_perdew86.o xc_perdew91.o \
	xc_dirac.o	xc_vwn.o	xc_lyp88.o			\
	eumac.o	\
	xc_tabcd.o	setACmat.o xc_rhogen.o\
	scat_mat.o	 \
	qpts.o		gausquad.o			   \
	a_dist.o	qdist.o		 xc_delleyw.o xc_lebedev.o
else
      OBJ = xc_getv.o  xc_quadv1.o  qpts.o  xc_quadv0.o xc_quadv0_a.o xc_perdew91.o

      OBJ_OPTIMIZE = \
	a_dist.o	  xc_beckew.o \
	xc_perdew86.o	 xc_dirac.o	\
	xc_becke88.o 	xc_lyp88.o	xc_vwn.o \
	eumac.o	\
	xc_fitv.o			   xc_get2oi.o	\
	gausquad.o	xc_lebedev.o		   \
	qdist.o	\
	xc_rhogen.o	   	\
	setACmat.o		xc_tabcd.o	\
	gat_mat.o		   scat_mat.o xc_delleyw.o
endif

     LIBRARY = libnwdft.a
   LIB_DEFINES = -DDFT_TIMER

   USES_BLAS = xc_fitv.F xc_quadv0_a.F xc_quadv1.F xc_dirac.F xc_tabcd.F xc_rhogen.F

include ../../config/makefile.h
include ../../config/makelib.h


