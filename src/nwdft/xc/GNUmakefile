# $Id: GNUmakefile,v 1.35 1996-10-23 01:15:52 gutowski Exp $
ifeq ($(TARGET),SGITFP)
	OBJ = xc_getv.o xc_quadv0.o xc_quadv0_a.o xc_quadv1.o xc_fitv.o xc_get2oi.o qpts.o gat_mat.o	xc_beckew.o setACmat.o xc_neighb.o 
	OBJ_OPTIMIZE = xc_tabcd.o  xc_rhogen.o\
	xc_becke88.o	xc_perdew86.o xc_perdew91.o xc_pw91lda.o\
	xc_dirac.o	xc_vwn.o	xc_lyp88.o			\
	eumac.o	\
	scat_mat.o	 \
	qpts.o		gausquad.o			   \
	a_dist.o	qdist.o		 xc_delleyw.o xc_lebedev.o \
        xc_vwnrpa.o xc_vwn_1.o xc_vwn_2.o xc_vwn_3.o xc_vwn_4.o
else
      OBJ = xc_getv.o  xc_quadv1.o  qpts.o  xc_quadv0.o xc_quadv0_a.o xc_perdew91.o 


      OBJ_OPTIMIZE = xc_rhogen.o		xc_pw91lda.o   	\
	a_dist.o	  xc_beckew.o 	 xc_neighb.o\
	xc_perdew86.o	 xc_dirac.o	\
	xc_becke88.o 	xc_lyp88.o	xc_vwn.o \
	eumac.o	\
	xc_fitv.o			   xc_get2oi.o	\
	gausquad.o	xc_lebedev.o		   \
	qdist.o	\
	setACmat.o		xc_tabcd.o	\
	gat_mat.o		   scat_mat.o xc_delleyw.o \
        xc_vwnrpa.o xc_vwn_1.o xc_vwn_2.o xc_vwn_3.o xc_vwn_4.o
endif

     LIBRARY = libnwdft.a
 LIB_DEFINES = -DDFT_TIMER
LIB_INCLUDES = -I../include


   USES_BLAS = xc_fitv.F xc_quadv0_a.F xc_quadv1.F xc_dirac.F xc_tabcd.F xc_rhogen.F

include ../../config/makefile.h
include ../../config/makelib.h


# DO NOT EDIT BENEATH THIS LINE ... GENERATED AUTOMATICALLY
$(LIBRARY_PATH)(xc_tabcd.o):	$(INCDIR)/bas.fh 
$(LIBRARY_PATH)(xc_rhogen.o):	$(INCDIR)/bas.fh 
$(LIBRARY_PATH)(xc_quadv1.o):	$(INCDIR)/util.fh ../include/timer.fh $(INCDIR)/tcgmsg.fh $(INCDIR)/stdio.fh $(INCDIR)/rtdb.fh $(INCDIR)/msgids.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh 
$(LIBRARY_PATH)(xc_quadv0_a.o):	../include/timer.fh $(INCDIR)/tcgmsg.fh $(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(xc_quadv0.o):	$(INCDIR)/util.fh ../include/timer.fh $(INCDIR)/tcgmsg.fh $(INCDIR)/stdio.fh $(INCDIR)/rtdb.fh $(INCDIR)/msgids.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh 
$(LIBRARY_PATH)(xc_lebedev.o):	$(INCDIR)/tcgmsg.fh $(INCDIR)/stdio.fh $(INCDIR)/global.fh 
$(LIBRARY_PATH)(xc_getv.o):	$(INCDIR)/tcgmsg.fh $(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(xc_get2oi.o):	$(INCDIR)/tcgmsg.fh $(INCDIR)/global.fh 
$(LIBRARY_PATH)(xc_fitv.o):	../include/timer.fh $(INCDIR)/tcgmsg.fh $(INCDIR)/msgids.fh $(INCDIR)/global.fh 
$(LIBRARY_PATH)(xc_beckew.o):	$(INCDIR)/bas_staticP.fh 
