# $Id: GNUmakefile,v 1.1 1995-04-26 03:44:39 og845 Exp $
ifeq ($(NWCHEM_TARGET),SGITFP)
      OBJ_OPTIMIZE = dft_getvc.o
      OBJ = \
	dft_fitcd.o				   dft_fitvc.o	\
	dft_get2eri.o		
else
      OBJ = dft_fitcd.o
      OBJ_OPTIMIZE = \
	dft_fitvc.o	\
	dft_get2eri.o			   dft_getvc.o	
endif

      LIBRARY = libnwdft.a

      USES_BLAS = dft_fitcd.F dft_fitvc.F

include ../../config/makefile.h
include ../../config/makelib.h

