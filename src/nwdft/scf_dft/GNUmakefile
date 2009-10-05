# $Id: GNUmakefile,v 1.63 2007-10-12 00:05:41 niri Exp $
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
	diis_hamwgt.o \
	dft_fdist_init.o \
	cdft_util.o \
	dft_zora_utils.o 
	OBJ=\
	dft_swapab.o \
	dft_main0d.o \
	dft_scf.o \
	dft_fockbld.o \
	dft_scfutil.o \
	dft_symbr.o \
	dft_pstat.o\
        fukui.o \
	dft_canorg.o \
	dft_diagn.o \
	dft_dump_info.o \
	dft_fant.o

     LIBRARY = libnwdft.a

 LIB_DEFINES =

LIB_INCLUDES = -I../include -I../grid

   USES_BLAS = dft_mxovl.F diis_bld12.F diis_bsolve.F dft_densm.F dft_swapab.F dft_canorg.F

HEADERS = dftps.fh dft_fdist.fh

include ../../config/makefile.h
include ../../config/makelib.h
