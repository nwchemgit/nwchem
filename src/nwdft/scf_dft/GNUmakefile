# $Id: GNUmakefile,v 1.19 1996-01-23 00:29:01 d3e129 Exp $

      OBJ = mull_pop.o dft_scf.o dft_s2_value.o
      OBJ_OPTIMIZE = \
	diis_driver.o	diis_bsolve.o	\
	diis_core.o	diis_bld12.o diis_hamwgt.o diis_refvec.o	\
	dft_dencvg.o dft_scfcvg.o	dft_densm.o	\
	dft_invdiag.o	\
	dft_main0d.o dft_mxovl.o

     LIBRARY = libnwdft.a
 LIB_DEFINES = -DDFT_TIMER
LIB_INCLUDES = -I../include

   USES_BLAS = diis_bld12.F diis_bsolve.F dft_mxovl.F

include ../../config/makefile.h
include ../../config/makelib.h



