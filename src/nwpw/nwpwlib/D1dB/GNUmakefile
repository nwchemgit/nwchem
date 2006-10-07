#
# $Id: GNUmakefile,v 1.2 2006-10-07 00:10:04 bylaska Exp $
#

 OBJ_OPTIMIZE = D1dB.o Dne.o Dnex.o
 OBJ = 

ifeq ($(USE_MPIF),y)
  OBJ_OPTIMIZE += D1dB-mpi.o
else
  OBJ_OPTIMIZE += D1dB-tcgmsg.o
endif



      LIBRARY = libnwpwlib.a

    USES_BLAS =

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

