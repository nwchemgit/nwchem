# $Id$
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
	cdft_util.o
	OBJ=\
	dft_swapab.o \
	dft_main0d.o \
	dft_scf.o \
	print_irrep_occp.o \
	dft_fockbld.o \
	dft_scfutil.o \
	dft_symbr.o \
	dft_pstat.o\
	fukui.o \
	dft_canorg.o \
	dft_diagn.o \
	dft_dump_info.o \
	dft_fant.o\
	dft_printinfo.o\
        dft_frac_mo.o

     LIBRARY = libnwdft.a
ifndef OLD_GA
ifdef USE_SCALAPACK
ifeq ($(shell cat ${NWCHEM_TOP}/src/tools/build/config.h | awk '/HAVE_ELPA/ {print "Y"}'),Y)
 LIB_DEFINES += -DUSE_ELPA
endif
ifeq ($(shell nm ${NWCHEM_TOP}/src/tools/build/global/src/scalapack.o|egrep ga_pzheevd |awk '/ga_pzheevd/ { print "Y"}'),Y)
 LIB_DEFINES += -DUSE_GAPZHEEVD
endif
_PDSYEVR = $(shell cat ${NWCHEM_TOP}/src/tools/build/config.h | awk ' /HAVE_PDSYEVR\ 1/ {print "Y"}')
ifeq ($(_PDSYEVR),Y) 
 LIB_DEFINES += -DUSE_PDSYEVR
endif
endif
endif

LIB_INCLUDES = -I../include -I../grid

   USES_BLAS = dft_mxovl.F diis_bld12.F diis_bsolve.F dft_densm.F dft_swapab.F dft_canorg.F cdft_util.F

HEADERS = dftps.fh dft_fdist.fh dft_conv.fh

include ../../config/makefile.h
include ../../config/makelib.h
