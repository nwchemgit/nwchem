#
# $Id: GNUmakefile,v 1.7 2002-04-24 21:59:17 bylaska Exp $
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
                sdminimize.o \
		nominimize.o
                


	  OBJ = 



      LIBRARY = libpspw.a

    USES_BLAS = cgsdv5.F cgsd_energy.F
                

  LIB_DEFINES = 
  LIB_INCLUDES = -I../../nwpwlib/include


include ../../../config/makefile.h
include ../../../config/makelib.h

