# $Id: GNUmakefile,v 1.2 1994-11-04 04:23:55 og845 Exp $

      OBJ_OPTIMIZE = \
	a_dist.o			   becke.o	\
	becke88.o			   denrow.o	\
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
