
# $Id$
#
 OBJ_OPTIMIZE = \
		cpsd_input.o \
		cpmd_input.o \
                cpsdv5.o \
                cpmdv5.o \
		mmsdv1.o \
                inner_loop_md.o \
                mm_loop.o \
                psi_lmbda.o psi_lmbda_sic.o psi_lmbda_paw.o psi_lmbda2.o psi_lmbda3.o \
                psi_lmbda_omp.o \
		fei_output.o \
		dipole_motion.o \
		pspw_et.o \
		runsocket.o
##                cpmd_qmmm_start.o cpmd_qmmm_step.o cpmd_qmmm_stop.o inner_loop_qmmm_step.o


	  OBJ_OPTIMIZE += psi_lmbda_omp.o inner_loop.o


      LIBRARY = libpspw.a

    USES_BLAS = psi_lmbda.f psi_lmbda2.f psi_lmbda3.f psi_lmbda_sic.f psi_lmbda_paw.f \
                inner_loop.F inner_loop_md.F cpsdv5.F cpmdv5.F psi_lmbda2.f psi_lmbda_omp.F \
		mmsdv1.F inner_loop_qmmm_step.F cpmd_qmmm_start.F cpmd_qmmm_stop.F mm_loop.F pspw_et.F

  LIB_DEFINES =


include ../../../config/makefile.h
include ../../../config/makelib.h

