# $Id: GNUmakefile,v 1.5 1995-06-23 00:46:29 og845 Exp $
ifeq ($(NWCHEM_TARGET),SGITFP)
      OBJ_OPTIMIZE = dft_getvc.o
      OBJ = \
	dft_fitcd.o				   dft_fitvc.o	\
	dft_get2eri.o		dft_genspm.o dft_getspm.o
else
      OBJ = dft_fitcd.o dft_genspm.o dft_getspm.o
      OBJ_OPTIMIZE = \
	dft_fitvc.o	\
	dft_get2eri.o			   dft_getvc.o	
endif

      LIBRARY = libnwdft.a
   LIB_DEFINES = -DDFT_TIMER 


      USES_BLAS = dft_fitcd.F dft_fitvc.F dft_genspm.F

include ../../config/makefile.h
include ../../config/makelib.h

