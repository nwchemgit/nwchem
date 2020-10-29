#
# $Id$
#
 OBJ_OPTIMIZE = \
		cgsd_input.o \
		pspw_noit_energy.o \
		pspw_energy.o \
		pspw_gradient.o \
		pspw_bq_gradient.o \
		pspw_stress.o \
		cgsdv5.o \
		cgmontecarlov1.o \
		cgoptimize1.o \
		cgmdv5.o \
		mmmdv1.o \
		cgsd_noit_energy.o \
		cgsd_energy.o \
		cgsd_excited.o \
		cgsd_oep_eigenvalues.o \
                cgminimize.o \
		cgminimize2.o \
                bfgsminimize.o  \
                bfgsminimize2.o  \
                bfgsminimize3.o \
                bybminimize0.o \
                bybminimize.o \
                bybminimize2.o \
                sdminimize.o \
		nominimize.o \
		md_energy.o \
		kp_structure.o

	  OBJ = 


      LIBRARY = libpspw.a

    USES_BLAS = cgsdv5.F cgmdv5.F cgmontecarlov1.F cgoptimize1.F cgsd_energy.F \
            bybminimize.F bybminimize2.F pspw_bq_gradient.F \
            cgminimize.F mmmdv1.F bfgsminimize2.F bfgsminimize3.F md_energy.F \
            cgminimize2.F cgsd_noit_energy.F bfgsminimize.F 

  LIB_DEFINES = 
  LIB_INCLUDES = -I../../nwpwlib/include


include ../../../config/makefile.h
include ../../../config/makelib.h

