#
# $Id: GNUmakefile,v 1.46 1999-11-05 20:17:10 edo Exp $
#

   OBJ = dft_grad_top.o  dftg_force.o \
         dftg_getxc.o dft_gradients.o \
         grad_dft.o      dftg_quadv0.o     \
         dftg_getcden.o  dftg_getxden.o  \
         dftg_bal_sum.o
   OBJ_OPTIMIZE = dftg_cdfit.o  dftg_quadv0a.o  dftg_quadv0b.o  dftg_ssw.o\
                  dftg_grad2.o dftg_hf_x.o dint_block_j.o j_nucder_gen.o  \
                  dftg_beckew.o dftg_wderiv.o
ifeq ($(NWCHEM_TARGET),SGITFP)
   OBJ = dftg_cdfit.o  dft_grad_top.o  dftg_force.o \
         dftg_getxc.o dft_gradients.o dftg_grad2.o \
         grad_dft.o      dftg_quadv0.o     \
         dftg_getcden.o  dftg_getxden.o\
         dftg_bal_sum.o dftg_hf_x.o dftg_quadv0a.o dftg_beckew.o dftg_wderiv.o\
	 dint_block_j.o j_nucder_gen.o  dftg_ssw.o
   OBJ_OPTIMIZE = dftg_quadv0b.o 
endif

   USES_BLAS = dftg_cdfit.F dftg_quadv0a.F dftg_quadv0b.F dftg_wderiv.F j_nucder_gen.F dft_gradients.F dftg_force.F

LIB_INCLUDES = -I../nwdft/include

   LIBRARY = libdftgrad.a

include ../config/makefile.h
include ../config/makelib.h

