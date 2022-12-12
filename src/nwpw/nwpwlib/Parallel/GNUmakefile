#
# $Id$
#


 OBJ_OPTIMIZE = Parallel.o
 	OBJ = 
	HEADERS = stupid_mpi4.fh

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
  OBJ_OPTIMIZE += Parallel-mpi.o
else
  OBJ_OPTIMIZE += Parallel-tcgmsg.o
endif


      LIBRARY = libnwpwlib.a

    USES_BLAS = Parallel-mpi.F Parallel.F

  LIB_DEFINES = 


include ../../../config/makefile.h
include ../../../config/makelib.h


ifeq ($(USE_MPIF),y)
  LIB_DEFINES += -DMPI
  ifeq ($(DEV_GA),y)
    LIB_DEFINES += -DGA_DEV
  endif
else
  LIB_DEFINES += -DTCGMSG
endif

ifeq ($(USE_MPIF4),y)
  LIB_DEFINES += -DMPI4
endif

