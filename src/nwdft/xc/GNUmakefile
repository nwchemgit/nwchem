# $Id: GNUmakefile,v 1.11 1995-02-17 06:53:04 d3e129 Exp $
ifeq ($(NWCHEM_TARGET),SGITFP)
	OBJ = getvxc.o quadvxc0.o quadvxc1.o fitvxc.o gaget2oi.o perdew86.o
	OBJ_OPTIMIZE = \
	becke.o	\
	becke88.o			   \
	dirac.o				   eumac.o	\
	lyp88.o			\
	tabcd.o	setACmat.o rho_gen.o\
	vwn.o				   scat_mat.o \
	qpts.o		gausquad.o			   \
	a_dist.o	qdist.o		getvxc.o 
else
      OBJ = getvxc.o perdew86.o
      OBJ_OPTIMIZE = \
	a_dist.o			   becke.o	\
	becke88.o			   \
	dirac.o				   eumac.o	\
	fitvxc.o			   gaget2oi.o	\
	gausquad.o			   \
	lyp88.o				   qdist.o	\
	qpts.o				   quadvxc0.o	\
	quadvxc1.o			   rho_gen.o	\
	setACmat.o			   tabcd.o	\
	vwn.o				   scat_mat.o 
endif

     LIBRARY = libnwdft.a
   LIB_DEFINES = -DDFT_TIMER

   USES_BLAS = fitvxc.F quadvxc0.F quadvxc1.F dirac.F tabcd.F rho_gen.F

include ../../config/makefile.h
include ../../config/makelib.h


