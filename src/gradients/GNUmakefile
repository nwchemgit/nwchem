# $Id: GNUmakefile,v 1.28 1997-03-25 21:21:56 d3e129 Exp $

 OBJ = gradients.o grad_force.o grad1.o scf_gradient.o grad_getdens.o \
       grad_dens.o grad_store.o grad_inp.o
 OBJ_OPTIMIZE = 

ifeq ($(TARGET),SGITFP)
 OBJ += grad2.o
else
# OBJ_OPTIMIZE += grad2.
 OBJ += grad2.o
endif

   USES_BLAS = grad2.F

   LIBRARY = libgradients.a

include ../config/makefile.h
include ../config/makelib.h

