#
# $Id: GNUmakefile,v 1.14 2004-03-09 21:25:57 bylaska Exp $
#
 OBJ_OPTIMIZE = \
		cgsd_input.o \
		pspw_noit_energy.o \
		pspw_energy.o \
		pspw_gradient.o \
		pspw_stress.o \
		cgsdv5.o \
		cgsd_noit_energy.o \
		cgsd_energy.o \
		cgsd_excited.o \
		cgsd_oep_eigenvalues.o \
                cgminimize.o \
		cgminimize2.o \
                bfgsminimize.o  \
                bfgsminimize2.o  \
                bfgsminimize3.o \
                bybminimize.o \
                sdminimize.o \
		nominimize.o


	  OBJ = 



      LIBRARY = libpspw.a

    USES_BLAS = cgsdv5.F cgsd_energy.F

  LIB_DEFINES = 
  LIB_INCLUDES = -I../../nwpwlib/include


include ../../../config/makefile.h
include ../../../config/makelib.h

