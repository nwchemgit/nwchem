#
# $Id: GNUmakefile,v 1.11 2003-04-22 18:05:10 bylaska Exp $
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
                cgminimize.o \
		cgminimize2.o \
                bfgsminimize.o \
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

