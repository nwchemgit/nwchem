# $Id: GNUmakefile,v 1.15 2004-08-26 23:35:11 edo Exp $

SUBDIRS = dft 

OBJ_OPTIMIZE =  hess_anal.o \
                hess_init.o \
                nucdd_cont.o \
                nucdd_calc.o \
                oned_calc.o \
                onedd_cont.o \
                get_wdens.o \
                onedd_calc.o \
                hess_dipder.o \
                twodd_cont.o \
                twodd_coul_ex.o \
                shell_fock_build.o \
                hess_cphf.o \
                hess_pxfy.o \
                hess_sxsy.o \
                sx_dra.o \
                hess_wrt.o \
                hess_restart.o \
                hess_hssout.o hess_print.o \
                hess_tidy.o

     LIBRARY = libhessian.a

   HEADERS = hess_info.fh

LIB_INCLUDES = -I../../../nwdft/include

   USES_BLAS = 	 \
    onedd_calc.F    hess_restart.F shell_fock_build.F  twodd_cont.F

include ../../config/makefile.h
include ../../config/makelib.h

