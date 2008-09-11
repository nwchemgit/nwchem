#
# $Id: GNUmakefile,v 1.18 2008-09-11 21:26:50 bylaska Exp $
#

 OBJ_OPTIMIZE = \
		dipole.o \
		orthocheck.o \
		psi.o \
		psi_spin2.o \
		geodesic.o \
		geodesic2.o \
		Grassman.o \
		psi_read.o \
		psi_write.o \
		v_psi_read.o \
		v_psi_write.o \
		psi_new.o \
		v_psi_new.o \
		psi_history.o \
		psi_write_filename.o \
		psi_read_filename.o \
		berry_phase_pol.o




      LIBRARY = libpspw.a

    USES_BLAS = psi.F geodesic.F geodesic2.F Grassman.F dipole.F psi_history.F skew.F geodesic1.F berry_phase_pol.F

  LIB_DEFINES = 


include ../../../../config/makefile.h
include ../../../../config/makelib.h

ifeq ($(_G77V33),Y)
FOPTIMIZE += -O0
endif

