#
# $Id: GNUmakefile,v 1.10 1998-09-14 23:06:42 d3e129 Exp $
#

 SUBDIRS = analytic/dft

 OBJ = scf_hessian.o rhf_hessian.o hess_quad.o hess_blas.o
 
 OBJ_OPTIMIZE =  

   USES_BLAS = hess_blas.F 

   LIBRARY = libhessian.a

include ../config/makefile.h

include ../config/makelib.h

ifeq ($(TARGET),KSR)
hess_quad.o:    hess_quad.F
	f77 -c -g  -I.  -I/home4/d3e129/nwdev/src/include  -DKSR -DPARALLEL_DIAG -DEXT_INT  hess_quad.F
endif

