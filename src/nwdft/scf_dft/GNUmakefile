# $Id: GNUmakefile,v 1.25 1996-10-01 21:38:02 d3h449 Exp $

      OBJ = dft_scf.o dft_s2_value.o dft_main0d.o
      OBJ_OPTIMIZE = \
	diis_driver.o	diis_bsolve.o	\
	diis_core.o	diis_bld12.o diis_hamwgt.o diis_refvec.o	\
	dft_dencvg.o dft_scfcvg.o	dft_densm.o	\
	dft_invdiag.o	dft_init_guess.o \
	dft_mxovl.o

     LIBRARY = libnwdft.a
 LIB_DEFINES = -DDFT_TIMER
LIB_INCLUDES = -I../include

   USES_BLAS = diis_bld12.F diis_bsolve.F dft_mxovl.F dft_init_guess.F

include ../../config/makefile.h
include ../../config/makelib.h



# DO NOT EDIT BENEATH THIS LINE ... GENERATED AUTOMATICALLY
$(LIBRARY_PATH)(diis_refvec.o):	$(INCDIR)/tcgmsg.fh 
$(LIBRARY_PATH)(diis_hamwgt.o):	$(INCDIR)/mafdecls.fh 
$(LIBRARY_PATH)(diis_driver.o):	../include/timer.fh $(INCDIR)/tcgmsg.fh $(INCDIR)/global.fh 
$(LIBRARY_PATH)(diis_core.o):	$(INCDIR)/util.fh $(INCDIR)/tcgmsg.fh $(INCDIR)/stdio.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
$(LIBRARY_PATH)(diis_bsolve.o):	$(INCDIR)/tcgmsg.fh $(INCDIR)/stdio.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
$(LIBRARY_PATH)(diis_bld12.o):	$(INCDIR)/tcgmsg.fh $(INCDIR)/stdio.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
$(LIBRARY_PATH)(dft_scfcvg.o):	$(INCDIR)/util.fh $(INCDIR)/tcgmsg.fh $(INCDIR)/rtdb.fh $(INCDIR)/msgids.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
$(LIBRARY_PATH)(dft_scf.o):	$(INCDIR)/util.fh ../include/timer.fh $(INCDIR)/stdio.fh $(INCDIR)/rtdb.fh $(INCDIR)/msgids.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh ../include/cdft.fh 
$(LIBRARY_PATH)(dft_s2_value.o):	$(INCDIR)/tcgmsg.fh $(INCDIR)/stdio.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
$(LIBRARY_PATH)(dft_mxovl.o):	$(INCDIR)/tcgmsg.fh $(INCDIR)/msgids.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
$(LIBRARY_PATH)(dft_main0d.o):	$(INCDIR)/util.fh ../include/timer.fh $(INCDIR)/tcgmsg.fh $(INCDIR)/stdio.fh $(INCDIR)/schwarz.fh $(INCDIR)/rtdb.fh $(INCDIR)/msgids.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh ../include/cdft.fh 
$(LIBRARY_PATH)(dft_invdiag.o):	$(INCDIR)/tcgmsg.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
$(LIBRARY_PATH)(dft_dencvg.o):	$(INCDIR)/tcgmsg.fh $(INCDIR)/msgids.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
