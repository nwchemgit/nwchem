# $Id: GNUmakefile,v 1.8 1995-01-25 23:18:38 og845 Exp $

      OBJ = mull_pop.o dft.o 
      OBJ_OPTIMIZE = \
	bsolve.o	\
	diis.o		dencvg.o scfcvg.o gasbld.o\
	diisdriver.o			   gacond.o	\
	gafact.o	\
	hamwght.o	\
	refvec.o			   scfban.o	\
	scfdriver.o

     LIBRARY = libnwdft.a
   LIB_DEFINES = -DDFT_TIMER

   USES_BLAS = gasbld.F bsolve.F

include ../../config/makefile.h
include ../../config/makelib.h



