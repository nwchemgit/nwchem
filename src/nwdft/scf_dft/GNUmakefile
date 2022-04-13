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
	dft_fock_efield.o\
	dft_frac_mo.o  dft_nwmolden.o dft_spinflip.o dft_quickguess.o

     LIBRARY = libnwdft.a
#_USE_SCALAPACK = $(shell cat ${NWCHEM_TOP}/src/tools/build/config.h | awk ' /HAVE_SCALAPACK\ 1/ {print "Y"}')


LIB_INCLUDES = -I../include -I../grid

   USES_BLAS = dft_mxovl.F diis_bld12.F diis_bsolve.F dft_densm.F dft_swapab.F dft_canorg.F cdft_util.F dft_diagn.F dft_spinflip.F

HEADERS = dftps.fh dft_fdist.fh dft_conv.fh

include ../../config/makefile.h
ifeq ($(_USE_SCALAPACK),Y)
_USE_ELPA = $(shell ${GA_PATH}/bin/ga-config  --use_elpa| awk ' /1/ {print "Y"}')

ifeq ($(_USE_ELPA),Y)
 LIB_DEFINES += -DUSE_ELPA
endif
#assuming we have GAPZHEEVD
GOTGAPZHEEVD= Y
ifeq ($(GOTGAPZHEEVD),Y)
 LIB_DEFINES += -DUSE_GAPZHEEVD
endif
_PDSYEVR = $(shell [ -e ${NWCHEM_TOP}/src/tools/install/lib/libga.a ] &&  nm ${NWCHEM_TOP}/src/tools/install/lib/libga.a|egrep T |egrep  ga_pdsyevr |awk '/ga_pdsyevr/ {print "Y"}')
ifeq ($(_PDSYEVR),Y) 
 LIB_DEFINES += -DUSE_PDSYEVR
endif
endif
include ../../config/makelib.h
