
 OBJ = scf_hessian.o rhf_hessian.o hess_quad.o hess_blas.o
 
 OBJ_OPTIMIZE =  

   USES_BLAS = hess_blas.F 

   LIBRARY = libhessian.a

include ../config/makefile.h

ifeq ($(TARGET),KSR)
hess_quad.o:    hess_quad.F
	f77 -c -g  -I.  -I/home4/d3e129/nwdev/src/include  -DKSR -DPARALLEL_DIAG -DEXT_INT  hess_quad.F
endif
ifeq ($(TARGET),CRAY-T3E)
hess_quad.o:    hess_quad.F
	$(MAKE) hess_quad.f
	f90 -c -g -F -O scalar1 hess_quad.f
	rm -f hess_quad.f
endif

include ../config/makelib.h

