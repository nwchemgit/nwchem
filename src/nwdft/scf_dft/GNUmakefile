# $Id: GNUmakefile,v 1.7 1994-12-31 01:22:05 og845 Exp $

      OBJ = mull_pop.o dft.o 
      OBJ_OPTIMIZE = \
	bsolve.o	\
	diis.o		scfcvg.o gasbld.o\
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



