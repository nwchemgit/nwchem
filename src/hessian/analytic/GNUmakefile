# $Id: GNUmakefile,v 1.12 2001-04-05 22:25:06 windus Exp $

SUBDIRS = dft 

OBJ_OPTIMIZE =  hess_anal.o \
                hess_init.o \
                nucdd_cont.o \
                nucdd_calc.o \
                oned_calc.o \
                onedd_cont.o \
                get_wdens.o \
                onedd_calc.o \
                twodd_cont.o \
                twodd_coul_ex.o \
                shell_fock_build.o \
                hess_cphf.o \
                hess_pxfy.o \
                hess_sxsy.o \
                sx_dra.o \
                hess_wrt.o \
                hess_restart.o \
                hess_hssout.o \
                hess_tidy.o

     LIBRARY = libhessian.a

   HEADERS = hess_info.fh

LIB_INCLUDES = -I../../../nwdft/include

   USES_BLAS = 	xc_d2expl.F \
		xc_rhogen_nucder.F \
		xc_rhogen_pert.F \
		xc_tabcd_cpks_lhs.F \
		xc_tabcd_cpks_rhs.F \
    onedd_calc.F    hess_restart.F shell_fock_build.F  twodd_cont.F

include ../../config/makefile.h
include ../../config/makelib.h

