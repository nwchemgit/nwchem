# $Id: GNUmakefile,v 1.31 1999-06-17 17:02:33 edo Exp $

 OBJ = gradients.o grad_force.o grad1.o scf_gradient.o  \
       grad_dens.o grad_inp.o ga_reorder.o
 OBJ_OPTIMIZE = grad2.o grad_getdens.o

#ifeq ($(TARGET),SGITFP)
# OBJ += grad2.o
#else
# OBJ_OPTIMIZE += grad2.
#endif

   USES_BLAS = grad2.F ga_reorder.F grad_dens.F

   LIBRARY = libgradients.a

include ../config/makefile.h
include ../config/makelib.h

