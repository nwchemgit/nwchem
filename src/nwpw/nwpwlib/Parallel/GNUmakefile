#
# $Id: GNUmakefile,v 1.7 2007-04-12 16:52:09 bylaska Exp $
#


 OBJ_OPTIMIZE = Parallel.o
 	OBJ = 
	HEADERS = stupid_mpi4.fh

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

ifeq ($(USE_MPIF4),y)
  LIB_DEFINES += -DMPI4
endif

