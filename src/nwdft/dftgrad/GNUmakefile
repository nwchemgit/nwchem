#
# $Id$
#

   OBJ = dft_grad_top.o \
         dftg_getxc.o dft_gradients.o \
         grad_dft.o          \
         dftg_getcden.o  dftg_getxden.o dftg_gridpart.o\
         dftg_gridv0a.o dftg_gridv0.o dftg_gridv0b.o \
         dftg_griddo.o
   OBJ_OPTIMIZE = dftg_cdfit.o dftg_quadv0b.o   \
                  dftg_grad2.o dint_block_j.o j_nucder_gen.o  \
                  dftg_beckew.o dftg_wderiv.o  dftg_gridssw.o dftg_gridsswold.o dftg_signf.o \
		  dftg_quadv0b_adft.o

   USES_BLAS = dftg_cdfit.F dftg_quadv0b.F dftg_wderiv.F j_nucder_gen.F dft_gradients.F   dftg_gridv0a.F dftg_gridv0b.F  dftg_quadv0b_adft.F

LIB_INCLUDES = -I../include -I../grid

   LIBRARY = libnwdft.a

include ../../config/makefile.h
include ../../config/makelib.h

