# $Id: GNUmakefile,v 1.3 1994-11-11 22:24:38 og845 Exp $

      OBJ_OPTIMIZE = \
	a_dist.o			   becke.o	\
	becke88.o			   \
	dirac.o				   eumac.o	\
	fitvxc.o			   gaget2oi.o	\
	gausquad.o			   getvxc.o	\
	lyp88.o				   qdist.o	\
	qpts.o				   quadvxc0.o	\
	quadvxc1.o			   rho_gen.o	\
	setACmat.o			   tabcd.o	\
	vwn.o				   vx_exact.o	

     LIBRARY = libnwdft.a

include ../../config/makefile.h
include ../../config/makelib.h

ifeq ($(TARGET),SGITFP)
fitvxc.o:	fitvxc.F
	f77   $(FOPTIONS) $(INCLUDES) $(DEFINES)  -c fitvxc.F -o fitvxc.o
endif
