#
# $Id: GNUmakefile,v 1.6 2001-12-21 16:50:29 bylaska Exp $
#
 OBJ_OPTIMIZE = \
				cgsd_input.o \
        		pspw_energy.o \
 				pspw_gradient.o \
 				pspw_stress.o \
				cgsdv5.o \
				cgsd_energy.o \
                cgminimize.o \
                bfgsminimize.o \
                sdminimize.o 
                


	  OBJ = 



      LIBRARY = libpspw.a

    USES_BLAS = cgsdv5.F cgsd_energy.F
                

  LIB_DEFINES = 
  LIB_INCLUDES = -I../../nwpwlib/include


include ../../../config/makefile.h
include ../../../config/makelib.h

