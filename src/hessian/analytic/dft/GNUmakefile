# $Id: GNUmakefile,v 1.2 1998-09-16 08:09:54 d3e129 Exp $

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

   USES_BLAS = 	xc_d2expl.F \
		xc_rhogen_nucder.F \
		xc_rhogen_pert.F \
		xc_tabcd_cpks_lhs.F \
		xc_tabcd_cpks_rhs.F

include ../../../config/makefile.h
include ../../../config/makelib.h

