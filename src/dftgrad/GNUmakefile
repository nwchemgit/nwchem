   OBJ = dft_grad_top.o  dftg_force.o \
         dftg_getxc.o       dftg_grad2.o \
         grad_dft.o      dftg_quadv0.o     \
         dftg_getcden.o  dftg_getxden.o  \
         dftg_bal_sum.o dftg_hf_x.o dftg_beckew.o dftg_wderiv.o
   OBJ_OPTIMIZE = dftg_cdfit.o  dftg_quadv0a.o  dftg_quadv0b.o  
ifeq ($(NWCHEM_TARGET),SGITFP)
   OBJ = dftg_cdfit.o  dft_grad_top.o  dftg_force.o \
         dftg_getxc.o       dftg_grad2.o \
         grad_dft.o      dftg_quadv0.o     \
         dftg_getcden.o  dftg_getxden.o\
         dftg_bal_sum.o dftg_hf_x.o dftg_quadv0a.o dftg_beckew.o dftg_wderiv.o
   OBJ_OPTIMIZE = dftg_quadv0b.o 
endif
#old   USES_BLAS = dftg_force.F dftg_quadv0a.F dftg_cdfit.F dftg_quadv0b.F dftg_wderiv.F dftg_wderiv.F
   USES_BLAS = dftg_cdfit.F dftg_quadv0a.F dftg_quadv0b.F dftg_wderiv.F
LIB_INCLUDES = -I../nwdft/include

   LIBRARY = libdftgrad.a

include ../config/makefile.h
include ../config/makelib.h

