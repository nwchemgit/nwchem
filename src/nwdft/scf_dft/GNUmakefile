# $Id: GNUmakefile,v 1.10 1995-05-10 03:52:45 og845 Exp $

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



