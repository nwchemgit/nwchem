# $Id: GNUmakefile,v 1.15 1996-05-01 23:37:05 og845 Exp $
ifeq ($(TARGET),SGITFP)
      OBJ_OPTIMIZE = dft_getvc.o dft_n3c.o
      OBJ =  dft_3cincor.o\
	dft_fitcd.o				   dft_fitvc.o	\
	dft_get2eri.o		dft_genspm.o dft_getspm.o 
else
      OBJ = dft_genspm.o dft_getspm.o 
      OBJ_OPTIMIZE = dft_fitcd.o dft_3cincor.o dft_n3c.o \
	dft_fitvc.o	\
	dft_get2eri.o			   dft_getvc.o	
endif

      LIBRARY = libnwdft.a
  LIB_DEFINES = -DDFT_TIMER 
 LIB_INCLUDES = -I../include



      USES_BLAS = dft_fitcd.F dft_fitvc.F dft_genspm.F

include ../../config/makefile.h
include ../../config/makelib.h

# DO NOT EDIT BENEATH THIS LINE ... GENERATED AUTOMATICALLY
$(LIBRARY_PATH)(dft_3cincor.o):	$(INCDIR)/tcgmsg.fh $(INCDIR)/schwarz.fh $(INCDIR)/rtdb.fh $(INCDIR)/msgids.fh $(INCDIR)/global.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(dft_fitcd.o):	$(INCDIR)/util.fh ../include/timer.fh $(INCDIR)/tcgmsg.fh $(INCDIR)/schwarz.fh $(INCDIR)/rtdb.fh $(INCDIR)/msgids.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/context.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(dft_getvc.o):	$(INCDIR)/tcgmsg.fh $(INCDIR)/stdio.fh $(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(dft_vc_solve_init.o):	$(INCDIR)/tcgmsg.fh 
$(LIBRARY_PATH)(dft_genspm.o):	$(INCDIR)/tcgmsg.fh $(INCDIR)/stdio.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(dft_fitvc.o):	../include/timer.fh $(INCDIR)/tcgmsg.fh $(INCDIR)/schwarz.fh $(INCDIR)/rtdb.fh $(INCDIR)/msgids.fh $(INCDIR)/global.fh $(INCDIR)/context.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(dft_n3c.o):	$(INCDIR)/tcgmsg.fh $(INCDIR)/rtdb.fh $(INCDIR)/msgids.fh $(INCDIR)/global.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(dft_vc_solve.o):	$(INCDIR)/tcgmsg.fh $(INCDIR)/schwarz.fh $(INCDIR)/rtdb.fh $(INCDIR)/msgids.fh $(INCDIR)/global.fh 
$(LIBRARY_PATH)(dft_getspm.o):	$(INCDIR)/tcgmsg.fh $(INCDIR)/stdio.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh 
$(LIBRARY_PATH)(dft_get2eri.o):	$(INCDIR)/tcgmsg.fh $(INCDIR)/global.fh $(INCDIR)/bas.fh 
