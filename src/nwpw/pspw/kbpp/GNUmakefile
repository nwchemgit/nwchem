#
# $Id: GNUmakefile,v 1.13 2008-03-17 17:53:57 bylaska Exp $
#
 OBJ_OPTIMIZE = \
                psp_formatter_input.o \
                psp_formatter.o \
                psp_stress_formatter.o \
		pawppv1.o \
	      	hghppv1.o \
		HGH_local.o \
		HGH_nonlocal.o \
		Tesseral.o \
	      	kbppv3.o \
	      	kbppv3-stress.o \
                setup_kbpp.o \
		kbpp_ray.o \
                integrate_kbppv3.o \
                integrate_kbppv3_new.o \
                integrate_kbppv3_ray.o \
                integrate_stress.o \
                integrate_stress_new.o \
                integrate_stress_ray.o \
                aorbs_formatter.o \
	      		kbpp_orb.o \
                integrate_kbpp_orb.o  \
                integrate_pawppv1.o \
                integrate_paw_ray.o \
		kbppv3d.o kbppv3d_stress.o \
		integrate_kbppv3d.o integrate_d_stress.o \
		integrate_kbppv3d_new.o integrate_d_stress_new.o

	  OBJ = 

      LIBRARY = libpspw.a

    USES_BLAS = setup_kbpp.f integrate_kbppv3.f integrate_kbppv3d.f integrate_kbpp_orb.F integrate_stress.F integrate_d_stress.F HGH_local.F HGH_nonlocal.F Tesseral.F hghppv1.F  pawppv1.F integrate_kbppv3d_new.f integrate_kbppv3_new.f integrate_d_stress_new.F

  LIB_DEFINES = 


include ../../../config/makefile.h
include ../../../config/makelib.h

