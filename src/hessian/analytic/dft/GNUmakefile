# $Id$

OBJ =  		xc_d2expl.o 
OBJ_OPTIMIZE =  xc_cpks_coeff.o \
		xc_exit_index.o \
		xc_init_index.o \
		xc_nucder_gen.o \
		dft_3dacc.o 

     LIBRARY = libhessian.a

# LIB_DEFINES = -DDFT_TIMER
 LIB_DEFINES =

LIB_INCLUDES = -I../../../nwdft/include

   USES_BLAS = 	xc_cpks_coeff.F \
                xc_d2expl.F \
                xc_nucder_gen.F

include ../../../config/makefile.h
include ../../../config/makelib.h

      FOPTIONS+= $(FOPTIMIZE)
      FDEBUG=
  ifeq ($(_FC),ifc) 
    ifeq ($(_IFCV8),Y)
      FOPTIONS+= -ansi_alias-
    endif	
  endif
  ifeq ($(_FC),efc) 
    ifeq ($(_IFCV8),Y)
      FOPTIONS+= -ansi_alias-
    endif	
  endif


