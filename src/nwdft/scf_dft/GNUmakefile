# $Id: GNUmakefile,v 1.47 2000-08-02 00:52:23 edo Exp $
#
#      OBJ = dft_scf.o dft_s2_value.o dft_main0d.o
#      OBJ_OPTIMIZE = \
#	diis_driver.o	diis_bsolve.o	\
#	diis_core.o	diis_bld12.o diis_hamwgt.o 	\
#	dft_dencvg.o dft_scfcvg.o	dft_densm.o	\
#	dft_invdiag.o	dft_init_guess.o \
#	dft_mxovl.o  dft_mxspin_ovlp.o
#

	OBJ_OPTIMIZE = \
	dft_dencvg.o \
	dft_densm.o \
	dft_invdiag.o \
	dft_mxovl.o \
	dft_mxspin_ovlp.o \
	dft_s2_value.o \
	dft_scfcvg.o \
	diis_bld12.o \
	diis_bsolve.o \
	diis_core.o \
	diis_driver.o \
	diis_hamwgt.o 
	OBJ=\
	dft_swapab.o \
	dft_main0d.o \
	dft_scf.o \
	dft_fockbld.o \
	dft_diag.o \
	dft_scfutil.o \
        fukui.o


     LIBRARY = libnwdft.a

# LIB_DEFINES = -DDFT_TIMER
 LIB_DEFINES =

LIB_INCLUDES = -I../include -I../grid

   USES_BLAS = dft_mxovl.F diis_bld12.F diis_bsolve.F dft_densm.F dft_swapab.F 

include ../../config/makefile.h
include ../../config/makelib.h


