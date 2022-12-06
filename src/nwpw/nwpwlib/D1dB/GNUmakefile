#
# $Id$
#

 OBJ_OPTIMIZE = D1dB.o Dne.o Dnex.o skew.o 
 OBJ_OPTIMIZE += Dne_omp.o

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
  OBJ_OPTIMIZE += D1dB-mpi.o D1dBs-mpi.o
else
  OBJ_OPTIMIZE += D1dB-tcgmsg.o
endif



      LIBRARY = libnwpwlib.a

    USES_BLAS = D1dB-mpi.F D1dBs-mpi.F Dne.F skew.F Dne_omp.F Dnex.F

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

