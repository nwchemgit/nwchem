# $Id: GNUmakefile,v 1.1 1998-09-09 13:21:12 bjohnson Exp $

OBJ_OPTIMIZE =  xc_cpks_coeff.o \
		xc_d2expl.o \
		xc_exit_index.o \
		xc_init_index.o \
		xc_nucder_gen.o \
		xc_rhogen_nucder.o \
		xc_rhogen_pert.o \
		xc_tabcd_cpks_lhs.o \
		xc_tabcd_cpks_rhs.o

     LIBRARY = libhessian.a

# LIB_DEFINES = -DDFT_TIMER
 LIB_DEFINES =

LIB_INCLUDES = -I../../../nwdft/include

   USES_BLAS = 	xc_d2expl.o \
		xc_rhogen_nucder.o \
		xc_rhogen_pert.o \
		xc_tabcd_cpks_lhs.o \
		xc_tabcd_cpks_rhs.o

include ../../../config/makefile.h
include ../../../config/makelib.h

