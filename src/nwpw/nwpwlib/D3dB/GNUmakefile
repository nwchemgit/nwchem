#
# $Id: GNUmakefile,v 1.5 2004-03-14 14:36:26 bylaska Exp $
#

 OBJ_OPTIMIZE = D3dB-new.o
 OBJ = 


      LIBRARY = libnwpwlib.a

    USES_BLAS = D3dB-new.F

  LIB_DEFINES = 



include ../../../config/makefile.h
ifeq ($(FC),pgf90)
 OBJ_OPTIMIZE = 
 OBJ = D3dB-new.o
endif
include ../../../config/makelib.h


ifeq ($(USE_MLIB),y)
  FOPTIONS += -DMLIB
endif

ifeq ($(USE_MPI),y)
  FOPTIONS += -DMPI
  LIB_INCLUDES += -I$(MPI_INCLUDE)
else
  FOPTIONS += -DTCGMSG
endif
