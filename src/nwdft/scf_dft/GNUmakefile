# $Id: GNUmakefile,v 1.32 1997-02-17 20:52:21 d3g681 Exp $

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


