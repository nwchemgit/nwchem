#
# $Id: GNUmakefile,v 1.1 2001-08-30 00:25:11 edo Exp $
#
# intf and aux break on IBM SP xlhpf with optimization ... why? 
 OBJ_OPTIMIZE = \
                cgsd_input.o \
        		pspw_energy.o \
 				pspw_gradient.o \
 				pspw_stress.o \
				cgsdv5.o \
				cgsd_energy.o \
                cgminimize.o \
                linesearch.o \
				Order_Values.o
                


	  OBJ = 



      LIBRARY = libpspw.a

    USES_BLAS = cgsdv5.F cgsd_energy.F
                

  LIB_DEFINES =


include ../../../config/makefile.h
include ../../../config/makelib.h

