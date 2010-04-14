#
# $Id: GNUmakefile,v 1.4 2007-04-11 23:41:14 bylaska Exp $
#

  OBJ_OPTIMIZE = CMatrix.o
  OBJ = 

ifeq ($(USE_MPIF),y)
  OBJ_OPTIMIZE += CMatrix-mpi.o
else
  OBJ_OPTIMIZE += CMatrix-tcgmsg.o
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
  ifneq ($(MPI_INCLUDE),)
    LIB_INCLUDES += -I$(MPI_INCLUDE)
  endif
else
  LIB_DEFINES += -DTCGMSG
endif

ifeq ($(USE_MPIF4),y)
  LIB_DEFINES += -DMPI4
endif



