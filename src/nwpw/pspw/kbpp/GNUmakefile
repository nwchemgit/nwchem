#
# $Id: GNUmakefile,v 1.5 2002-02-11 19:47:33 edo Exp $
#
 OBJ_OPTIMIZE = \
                psp_formatter_input.o \
                psp_formatter.o \
                psp_stress_formatter.o \
	      		kbppv3.o \
	      		kbppv3-stress.o \
                setup_kbpp.o \
                integrate_kbppv3.o \
                integrate_stress.o \
                aorbs_formatter.o \
	      		kbpp_orb.o \
                integrate_kbpp_orb.o 


	  OBJ = 

      LIBRARY = libpspw.a

    USES_BLAS = setup_kbpp.f integrate_kbppv3.f integrate_kbpp_orb.F integrate_stress.f

  LIB_DEFINES = 


include ../../../config/makefile.h
include ../../../config/makelib.h

