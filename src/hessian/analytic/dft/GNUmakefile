# $Id: GNUmakefile,v 1.4 1998-12-11 16:09:10 d3j191 Exp $

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

   USES_BLAS = 	xc_cpks_coeff.F \
                xc_d2expl.F \
		xc_rhogen_nucder.F \
		xc_rhogen_pert.F \
		xc_tabcd_cpks_lhs.F \
		xc_tabcd_cpks_rhs.F\
                xc_nucder_gen.F

include ../../../config/makefile.h
include ../../../config/makelib.h

