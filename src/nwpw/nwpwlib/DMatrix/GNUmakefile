#
# $Id: GNUmakefile,v 1.4 2007-04-11 23:41:14 bylaska Exp $
#

  OBJ_OPTIMIZE = DMatrix.o
  OBJ = 

ifeq ($(USE_MPIF),y)
  OBJ_OPTIMIZE += DMatrix-mpi.o
else
  OBJ_OPTIMIZE += DMatrix-tcgmsg.o
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

ifeq ($(USE_MPIF4),y)
  LIB_DEFINES += -DMPI4
endif



