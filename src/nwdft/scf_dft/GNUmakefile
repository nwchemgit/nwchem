# $Id: GNUmakefile,v 1.11 1995-10-17 08:55:02 d3g681 Exp $

      OBJ = mull_pop.o dft_scf.o dft_s2_value.o
      OBJ_OPTIMIZE = \
	diis_driver.o	diis_bsolve.o	\
	diis_core.o	diis_builds12.o diis_hamwgt.o diis_refvec.o	\
	dft_dencvg.o dft_scfcvg.o	dft_densm.o	\
	dft_invdiag.o	\
	dft_main0d.o

     LIBRARY = libnwdft.a
   LIB_DEFINES = -DDFT_TIMER

   USES_BLAS = diis_builds12.F diis_bsolve.F

include ../../config/makefile.h
include ../../config/makelib.h



# DO NOT EDIT BENEATH THIS LINE ... GENERATED AUTOMATICALLY
dft_dencvg.o :	 $(INCDIR)/tcgmsg.fh $(INCDIR)/msgids.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
dft_invdiag.o :	 $(INCDIR)/tcgmsg.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
dft_main0d.o :	 $(INCDIR)/timer.fh $(INCDIR)/tcgmsg.fh $(INCDIR)/stdio.fh $(INCDIR)/schwarz.fh $(INCDIR)/rtdb.fh $(INCDIR)/msgids.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh $(INCDIR)/bas.fh 
dft_s2_value.o :	 $(INCDIR)/tcgmsg.fh $(INCDIR)/stdio.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
dft_scf.o :	 $(INCDIR)/timer.fh $(INCDIR)/tcgmsg.fh $(INCDIR)/stdio.fh $(INCDIR)/rtdb.fh $(INCDIR)/msgids.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh $(INCDIR)/context.fh $(INCDIR)/bas.fh 
dft_scfcvg.o :	 $(INCDIR)/tcgmsg.fh $(INCDIR)/stdio.fh $(INCDIR)/rtdb.fh $(INCDIR)/msgids.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
diis_bsolve.o :	 $(INCDIR)/tcgmsg.fh $(INCDIR)/stdio.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
diis_builds12.o :	 $(INCDIR)/tcgmsg.fh $(INCDIR)/stdio.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
diis_core.o :	 $(INCDIR)/tcgmsg.fh $(INCDIR)/stdio.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
diis_driver.o :	 $(INCDIR)/timer.fh $(INCDIR)/tcgmsg.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
diis_hamwgt.o :	 $(INCDIR)/mafdecls.fh 
diis_refvec.o :	 $(INCDIR)/tcgmsg.fh $(INCDIR)/global.fh 
mull_pop.o :	 $(INCDIR)/tcgmsg.fh $(INCDIR)/stdio.fh $(INCDIR)/msgids.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh $(INCDIR)/bas.fh 
