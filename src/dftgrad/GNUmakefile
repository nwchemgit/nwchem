   OBJ = dft_grad_top.o  dftg_force.o \
         dftg_getxc.o       dftg_grad2.o \
         grad_dft.o      dftg_quadv0.o     \
         dftg_get_cden.o  \
         dftg_bal_sum.o
   OBJ_OPTIMIZE = dftg_cdfit.o dftg_quadv0_b.o  
ifeq ($(NWCHEM_TARGET),SGITFP)
   OBJ = dftg_cdfit.o  dft_grad_top.o  dftg_force.o \
         dftg_getxc.o       dftg_grad2.o \
         grad_dft.o      dftg_quadv0.o     \
         dftg_get_cden.o  \
         dftg_bal_sum.o
   OBJ_OPTIMIZE = dftg_quadv0_b.o  
endif
   USES_BLAS = dftg_force.F dftg_quadv0.F dftg_cdfit.F dftg_quadv0_b.F

   LIBRARY = libdftgrad.a

include ../config/makefile.h
include ../config/makelib.h

