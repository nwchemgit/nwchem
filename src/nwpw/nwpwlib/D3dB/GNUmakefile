#
# $Id: GNUmakefile,v 1.4 2003-10-26 20:47:24 bylaska Exp $
#

 OBJ_OPTIMIZE = D3dB-new.o
 OBJ = 


      LIBRARY = libnwpwlib.a

    USES_BLAS = D3dB-new.F

  LIB_DEFINES = 


ifeq ($(USE_MLIB),y)
  FOPTIONS += -DMLIB
endif

include ../../../config/makefile.h
ifeq ($(FC),pgf90)
 OBJ_OPTIMIZE = 
 OBJ = D3dB-new.o
endif
include ../../../config/makelib.h

