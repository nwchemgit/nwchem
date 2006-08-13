#
# $Id: GNUmakefile,v 1.5 2006-08-13 01:03:25 bylaska Exp $
#


 OBJ_OPTIMIZE = Parallel.o
	  OBJ = 

ifeq ($(USE_MPIF),y)
  OBJ_OPTIMIZE += Parallel-mpi.o
else
  OBJ_OPTIMIZE += Parallel-tcgmsg.o
endif


      LIBRARY = libnwpwlib.a

    USES_BLAS = 

  LIB_DEFINES = 


include ../../../config/makefile.h
include ../../../config/makelib.h


ifeq ($(USE_MPIF),y)
  LIB_DEFINES += -DMPI
  LIB_INCLUDES += -I$(MPI_INCLUDE)
else
  LIB_DEFINES += -DTCGMSG
endif

