#
# $Id: GNUmakefile,v 1.1 2001-08-30 00:38:58 edo Exp $
#
# intf and aux break on IBM SP xlhpf with optimization ... why? 
 OBJ_OPTIMIZE = \
                psp_formatter.o \
                psp_stress_formatter.o \
                psp_formatter_input.o \
	      		kbppv3.o \
	      		kbppv3-stress.o \
                setup_kbpp.o \
                integrate_kbppv3.o \
                integrate_stress.o


	  OBJ = 

      LIBRARY = libpspw.a

    USES_BLAS = setup_kbpp.f integrate_kbppv3.f

  LIB_DEFINES = 


include ../../../config/makefile.h
include ../../../config/makelib.h

