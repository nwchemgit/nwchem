# $Id: GNUmakefile,v 1.30 1997-02-13 23:21:06 mdupuis Exp $

      OBJ = dft_scf.o dft_s2_value.o dft_main0d.o
      OBJ_OPTIMIZE = \
	diis_driver.o	diis_bsolve.o	\
	diis_core.o	diis_bld12.o diis_hamwgt.o diis_refvec.o	\
	dft_dencvg.o dft_scfcvg.o	dft_densm.o	\
	dft_invdiag.o	dft_init_guess.o \
	dft_mxovl.o  dft_mxspin_ovlp.o

     LIBRARY = libnwdft.a

# LIB_DEFINES = -DDFT_TIMER
 LIB_DEFINES =

LIB_INCLUDES = -I../include

   USES_BLAS = diis_bld12.F diis_bsolve.F dft_mxovl.F dft_mxspin_ovlp.F dft_init_guess.F

include ../../config/makefile.h
include ../../config/makelib.h


# DO NOT EDIT BENEATH THIS LINE ... GENERATED AUTOMATICALLY
$(LIBRARY_PATH)(dft_invdiag.o):	$(INCDIR)/tcgmsg.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
$(LIBRARY_PATH)(dft_dencvg.o):	$(INCDIR)/tcgmsg.fh $(INCDIR)/msgids.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
$(LIBRARY_PATH)(dft_mxspin_ovlp.o):	$(INCDIR)/tcgmsg.fh $(INCDIR)/stdio.fh $(INCDIR)/msgids.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh ../include/cdft.fh 
$(LIBRARY_PATH)(dft_init_guess.o):	$(INCDIR)/util.fh $(INCDIR)/stdio.fh $(INCDIR)/rtdb.fh $(INCDIR)/msgids.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh ../include/cdft.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(diis_bld12.o):	$(INCDIR)/tcgmsg.fh $(INCDIR)/stdio.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
$(LIBRARY_PATH)(dft_main0d.o):	$(INCDIR)/util.fh $(INCDIR)/tcgmsg.fh $(INCDIR)/stdio.fh $(INCDIR)/schwarz.fh $(INCDIR)/rtdb.fh $(INCDIR)/msgids.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh ../include/cdft.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(dft_mxovl.o):	$(INCDIR)/tcgmsg.fh $(INCDIR)/msgids.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh ../include/cdft.fh 
$(LIBRARY_PATH)(dft_scf.o):	$(INCDIR)/util.fh $(INCDIR)/stdio.fh $(INCDIR)/rtdb.fh $(INCDIR)/msgids.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh ../include/cdft.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(diis_refvec.o):	$(INCDIR)/tcgmsg.fh 
$(LIBRARY_PATH)(diis_hamwgt.o):	$(INCDIR)/mafdecls.fh 
$(LIBRARY_PATH)(diis_driver.o):	$(INCDIR)/tcgmsg.fh $(INCDIR)/stdio.fh $(INCDIR)/msgids.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
$(LIBRARY_PATH)(diis_core.o):	$(INCDIR)/util.fh $(INCDIR)/tcgmsg.fh $(INCDIR)/stdio.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
$(LIBRARY_PATH)(diis_bsolve.o):	$(INCDIR)/tcgmsg.fh $(INCDIR)/stdio.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
$(LIBRARY_PATH)(dft_scfcvg.o):	$(INCDIR)/util.fh $(INCDIR)/tcgmsg.fh $(INCDIR)/rtdb.fh $(INCDIR)/msgids.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
$(LIBRARY_PATH)(dft_s2_value.o):	$(INCDIR)/tcgmsg.fh $(INCDIR)/stdio.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
