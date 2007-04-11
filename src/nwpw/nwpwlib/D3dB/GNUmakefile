#
# $Id: GNUmakefile,v 1.12 2007-04-11 23:41:14 bylaska Exp $
#

 OBJ_OPTIMIZE = D3dB-new.o generate_map_indexes.o hilbert.o olist.o trans.o D3dB_pfft.o
 OBJ = 

ifeq ($(USE_MPIF),y)
  OBJ_OPTIMIZE += D3dB-mpi.o
else
  OBJ_OPTIMIZE += D3dB-tcgmsg.o
endif

      LIBRARY = libnwpwlib.a

    USES_BLAS = D3dB-new.F D3dB_pfft.F D3dB-tcgmsg.F

  LIB_DEFINES = 


include ../../../config/makefile.h
include ../../../config/makelib.h


ifeq ($(USE_MLIB),y)
  LIB_DEFINES += -DMLIB
endif

ifeq ($(USE_MPIF),y)
  LIB_DEFINES += -DMPI
  LIB_INCLUDES += -I$(MPI_INCLUDE)
else
  LIB_DEFINES += -DTCGMSG
endif

ifeq ($(USE_MPIF4),y)
  LIB_DEFINES += -DMPI4
endif
