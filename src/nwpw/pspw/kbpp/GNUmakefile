#
# $Id: GNUmakefile,v 1.8 2007-08-27 15:31:33 bylaska Exp $
#
 OBJ_OPTIMIZE = \
                psp_formatter_input.o \
                psp_formatter.o \
                psp_stress_formatter.o \
	      		hghppv1.o \
			HGH_local.o \
			HGH_nonlocal.o \
			Tesseral.o \
	      		kbppv3.o \
	      		kbppv3-stress.o \
                setup_kbpp.o \
                integrate_kbppv3.o \
                integrate_kbppv3_new.o \
                integrate_kbppv3_ray.o \
                integrate_stress.o \
                aorbs_formatter.o \
	      		kbpp_orb.o \
                integrate_kbpp_orb.o 


	  OBJ = 

      LIBRARY = libpspw.a

    USES_BLAS = setup_kbpp.f integrate_kbppv3.f integrate_kbpp_orb.F integrate_stress.F HGH_local.F HGH_nonlocal.F Tesseral.F hghppv1.F  

  LIB_DEFINES = 


include ../../../config/makefile.h
include ../../../config/makelib.h

