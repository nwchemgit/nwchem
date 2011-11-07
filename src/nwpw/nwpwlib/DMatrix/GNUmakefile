#
# $Id$
#

  OBJ_OPTIMIZE = DMatrix.o
  OBJ = 

ifeq ($(USE_MPIF),y)
  OBJ_OPTIMIZE += DMatrix-mpi.o
else
  OBJ_OPTIMIZE += DMatrix-tcgmsg.o
endif



      LIBRARY = libnwpwlib.a

    USES_BLAS = DMatrix-mpi.F DMatrix.F DMatrix-tcgmsg.F

  LIB_DEFINES = 


include ../../../config/makefile.h
include ../../../config/makelib.h


ifeq ($(USE_MLIB),y)
  LIB_DEFINES += -DMLIB
endif

ifeq ($(USE_MPIF),y)
  LIB_DEFINES += -DMPI
  ifneq ($(MPI_INCLUDE),)
    LIB_INCLUDES += $(patsubst -I-I%,-I%,-I$(MPI_INCLUDE))
  endif
else
  LIB_DEFINES += -DTCGMSG
endif

ifeq ($(USE_MPIF4),y)
  LIB_DEFINES += -DMPI4
endif



