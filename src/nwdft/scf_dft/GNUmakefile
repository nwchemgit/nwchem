# $Id: GNUmakefile,v 1.17 1995-12-13 03:08:58 og845 Exp $

      OBJ = mull_pop.o dft_scf.o dft_s2_value.o
      OBJ_OPTIMIZE = \
	diis_driver.o	diis_bsolve.o	\
	diis_core.o	diis_bld12.o diis_hamwgt.o diis_refvec.o	\
	dft_dencvg.o dft_scfcvg.o	dft_densm.o	\
	dft_invdiag.o	\
	dft_main0d.o dft_mxovl.o

     LIBRARY = libnwdft.a
 LIB_DEFINES = -DDFT_TIMER
LIB_INCLUDES = -I../include

   USES_BLAS = diis_bld12.F diis_bsolve.F dft_mxovl.F

include ../../config/makefile.h
include ../../config/makelib.h



# DO NOT EDIT BENEATH THIS LINE ... GENERATED AUTOMATICALLY
$(LIBRARY_PATH)(dft_dencvg.o):	 $(INCDIR)/tcgmsg.fh $(INCDIR)/msgids.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
$(LIBRARY_PATH)(dft_invdiag.o):	 $(INCDIR)/tcgmsg.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
$(LIBRARY_PATH)(dft_main0d.o):	 ../include/timer.fh $(INCDIR)/tcgmsg.fh $(INCDIR)/stdio.fh $(INCDIR)/schwarz.fh $(INCDIR)/rtdb.fh $(INCDIR)/msgids.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh ../include/cdft.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(dft_max_ovl.o):	 $(INCDIR)/tcgmsg.fh $(INCDIR)/msgids.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(dft_s2_value.o):	 $(INCDIR)/tcgmsg.fh $(INCDIR)/stdio.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
$(LIBRARY_PATH)(dft_scf.o):	 ../include/timer.fh $(INCDIR)/tcgmsg.fh $(INCDIR)/stdio.fh $(INCDIR)/rtdb.fh $(INCDIR)/msgids.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh $(INCDIR)/context.fh $(INCDIR)/cfock.fh ../include/cdft.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(dft_scfcvg.o):	 $(INCDIR)/tcgmsg.fh $(INCDIR)/stdio.fh $(INCDIR)/rtdb.fh $(INCDIR)/msgids.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
$(LIBRARY_PATH)(diis_bld12.o):	 $(INCDIR)/tcgmsg.fh $(INCDIR)/stdio.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
$(LIBRARY_PATH)(diis_bsolve.o):	 $(INCDIR)/tcgmsg.fh $(INCDIR)/stdio.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
$(LIBRARY_PATH)(diis_builds12.o):	 $(INCDIR)/tcgmsg.fh $(INCDIR)/stdio.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
$(LIBRARY_PATH)(diis_core.o):	 $(INCDIR)/tcgmsg.fh $(INCDIR)/stdio.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
$(LIBRARY_PATH)(diis_driver.o):	 ../include/timer.fh $(INCDIR)/tcgmsg.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
$(LIBRARY_PATH)(diis_hamwgt.o):	 $(INCDIR)/mafdecls.fh 
$(LIBRARY_PATH)(diis_refvec.o):	 $(INCDIR)/tcgmsg.fh $(INCDIR)/global.fh 
$(LIBRARY_PATH)(mull_pop.o):	 $(INCDIR)/tcgmsg.fh $(INCDIR)/stdio.fh $(INCDIR)/msgids.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh $(INCDIR)/bas.fh 
