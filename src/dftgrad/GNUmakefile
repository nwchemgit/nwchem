   OBJ = dft_grad_top.o  force_dft.o \
         getdexc.o       grad1_dft.o     grad2_dft.o \
         grad_dft.o      quaddexc0.o     \
         get_coul_den.o  grad_e_xc.o \
         setACtrans.o    grad_e_xc_gca.o
   OBJ_OPTIMIZE = grad_cdfit.o grad_xc.o
   USES_BLAS = force_dft.F grad_e_xc.F quaddexc0.F grad_cdfit.F grad_xc.F

   LIBRARY = libdftgrad.a

include ../config/makefile.h
include ../config/makelib.h

