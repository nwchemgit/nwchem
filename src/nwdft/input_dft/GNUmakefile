# $Id: GNUmakefile,v 1.26 1997-04-25 23:11:43 d3h449 Exp $

      OBJ = \
	dft_inpana.o \
	dft_rdinput.o \
        dft_input.o \
        tol_inp.o \
        grid_inp.o \
        xc_inp.o \
        conv_inp.o \
	xc_setquad.o

     LIBRARY = libnwdft.a

 LIB_DEFINES = -DDFT_TIMER
# LIB_DEFINES =

LIB_INCLUDES = -I../include

   USES_BLAS = xc_setquad.F

include ../../config/makefile.h
include ../../config/makelib.h


