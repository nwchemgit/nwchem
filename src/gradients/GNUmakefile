# $Id: GNUmakefile,v 1.29 1997-06-24 15:56:24 d3g681 Exp $

 OBJ = gradients.o grad_force.o grad1.o scf_gradient.o  \
       grad_dens.o grad_inp.o ga_reorder.o
 OBJ_OPTIMIZE = grad2.o grad_getdens.o

#ifeq ($(TARGET),SGITFP)
# OBJ += grad2.o
#else
# OBJ_OPTIMIZE += grad2.
#endif

   USES_BLAS = grad2.F

   LIBRARY = libgradients.a

include ../config/makefile.h
include ../config/makelib.h

