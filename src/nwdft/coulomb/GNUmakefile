# $Id: GNUmakefile,v 1.3 1995-05-08 19:46:21 og845 Exp $
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

      USES_BLAS = dft_fitcd.F dft_fitvc.F dft_genspm.o

include ../../config/makefile.h
include ../../config/makelib.h

