
# $Id$
#
 OBJ_OPTIMIZE = \
		cpsd_input.o \
		cpmd_input.o \
                cpsdv5.o \
                cpmdv5.o \
		mmsdv1.o \
                inner_loop.o \
                inner_loop_md.o \
                mm_loop.o \
                psi_lmbda.o psi_lmbda_sic.o psi_lmbda_paw.o psi_lmbda2.o psi_lmbda3.o \
		fei_output.o \
		pspw_et.o
##                cpmd_qmmm_start.o cpmd_qmmm_step.o cpmd_qmmm_stop.o inner_loop_qmmm_step.o


	  OBJ = 



      LIBRARY = libpspw.a

    USES_BLAS = psi_lmbda.f psi_lmbda2.f psi_lmbda3.f psi_lmbda_sic.f psi_lmbda_paw.f \
                inner_loop.F inner_loop_md.F cpsdv5.F cpmdv5.F psi_lmbda2.f \
 		mmsdv1.F inner_loop_qmmm_step.F cpmd_qmmm_start.F cpmd_qmmm_stop.F mm_loop.F pspw_et.F

  LIB_DEFINES =


include ../../../config/makefile.h
include ../../../config/makelib.h

