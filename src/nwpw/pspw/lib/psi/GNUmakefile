#
# $Id: GNUmakefile,v 1.15 2004-05-14 17:43:30 edo Exp $
#

 OBJ_OPTIMIZE = \
		dipole.o \
		orthocheck.o \
		psi.o \
		psi_spin2.o \
		geodesic.o \
		geodesic2.o \
		skew.o \
		Grassman.o \
		psi_read.o \
		psi_write.o \
		v_psi_read.o \
		v_psi_write.o \
		psi_new.o \
		v_psi_new.o \
		psi_history.o \
		psi_write_filename.o \
		psi_read_filename.o




      LIBRARY = libpspw.a

    USES_BLAS = psi.F geodesic.F geodesic2.F Grassman.F dipole.F psi_history.F skew.F geodesic1.F

  LIB_DEFINES = 


include ../../../../config/makefile.h
include ../../../../config/makelib.h

ifeq ($(_G77V33),Y)
FOPTIMIZE += -O0
endif

