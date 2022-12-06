#
# $Id$
#

  OBJ_OPTIMIZE = DMatrix.o
  OBJ = 

### make defaults USE_MPIF and USE_MPIF4 ###
ifeq ($(USE_MPIF),n)
 USE_MPIF = n
else
 USE_MPIF = y
endif
ifeq ($(USE_MPIF4),n)
 USE_MPIF4 = n
else
 USE_MPIF4 = y
endif


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
else
  LIB_DEFINES += -DTCGMSG
endif

ifeq ($(USE_MPIF4),y)
  LIB_DEFINES += -DMPI4
endif



