#
# $Id: GNUmakefile,v 1.5 2001-12-17 23:31:15 bylaska Exp $
#
 OBJ_OPTIMIZE = \
				cgsd_input.o \
        		pspw_energy.o \
 				pspw_gradient.o \
 				pspw_stress.o \
				cgsdv5.o \
				cgsd_energy.o \
                cgminimize.o \
                sdminimize.o 
                


	  OBJ = 



      LIBRARY = libpspw.a

    USES_BLAS = cgsdv5.F cgsd_energy.F
                

  LIB_DEFINES = 
  LIB_INCLUDES = -I../../nwpwlib/include


include ../../../config/makefile.h
include ../../../config/makelib.h

