# $Id: GNUmakefile,v 1.7 1994-11-27 10:38:40 og845 Exp $
      OBJ = getvxc.o
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

     LIBRARY = libnwdft.a

   USES_BLAS = fitvxc.F quadvxc0.F quadvxc1.F dirac.F tabcd.F rho_gen.F

include ../../config/makefile.h
include ../../config/makelib.h

ifeq ($(TARGET),SGITFP)
quadvxc0.o:	quadvxc0.F
	f77   $(FOPTIONS)  $(INCLUDES) $(DEFINES)  -c quadvxc0.F -o quadvxc0.o
quadvxc1.o:	quadvxc1.F
	f77   $(FOPTIONS)  $(INCLUDES) $(DEFINES)  -c quadvxc1.F -o quadvxc1.o
fitvxc.o:	fitvxc.F
	f77   $(FOPTIONS) $(INCLUDES) $(DEFINES)  -c fitvxc.F -o fitvxc.o
gaget2oi.o:	gaget2oi.F
	f77   $(FOPTIONS) $(INCLUDES) $(DEFINES)  -c gaget2oi.F -o gaget2oi.o
endif
