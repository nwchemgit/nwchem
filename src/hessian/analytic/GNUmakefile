# $Id: GNUmakefile,v 1.5 2000-05-30 14:47:29 windus Exp $

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
                hess_cphf.o \
                hess_pxfy.o \
                sx_dra.o \
                hess_tidy.o

     LIBRARY = libhessian.a

   HEADERS = hess_info.fh

LIB_INCLUDES = -I../../../nwdft/include

   USES_BLAS = 	xc_d2expl.o \
		xc_rhogen_nucder.o \
		xc_rhogen_pert.o \
		xc_tabcd_cpks_lhs.o \
		xc_tabcd_cpks_rhs.o

include ../../config/makefile.h
include ../../config/makelib.h

