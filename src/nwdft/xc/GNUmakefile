# $Id: GNUmakefile,v 1.15 1995-04-09 08:09:26 og845 Exp $
ifeq ($(NWCHEM_TARGET),SGITFP)
	OBJ = getvxc.o quadvxc0.o quadvxc1.o fitvxc.o gaget2oi.o qpts.o gat_mat.o	becke.o
	OBJ_OPTIMIZE = \
	becke88.o	perdew86.o 		   \
	dirac.o				   eumac.o	\
	lyp88.o			\
	tabcd.o	setACmat.o rho_gen.o\
	vwn.o		   scat_mat.o	 \
	qpts.o		gausquad.o			   \
	a_dist.o	qdist.o		getvxc.o delley.o lebedev.o
else
      OBJ = getvxc.o  quadvxc1.o quadvxc0.o qpts.o

      OBJ_OPTIMIZE = \
	a_dist.o	 perdew86.o		   becke.o	\
	becke88.o 		    \
	dirac.o				   eumac.o	\
	fitvxc.o			   gaget2oi.o	\
	gausquad.o	lebedev.o		   \
	lyp88.o				   qdist.o	\
	rho_gen.o	   	\
	setACmat.o			   tabcd.o	\
	vwn.o		gat_mat.o		   scat_mat.o delley.o
endif

     LIBRARY = libnwdft.a
   LIB_DEFINES = -DDFT_TIMER

   USES_BLAS = fitvxc.F quadvxc0.F quadvxc1.F dirac.F tabcd.F rho_gen.F

include ../../config/makefile.h
include ../../config/makelib.h


