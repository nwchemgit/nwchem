#
# $Id: GNUmakefile,v 1.4 2006-03-16 18:38:27 bylaska Exp $
#

 OBJ_OPTIMIZE = C3dB-new.o c_trans.o
	  OBJ = 

#ifeq ($(USE_MPIF),y)
#  OBJ_OPTIMIZE += C3dB-mpi.o
#else
  OBJ_OPTIMIZE += C3dB-tcgmsg.o
#endif



      LIBRARY = libnwpwlib.a

    USES_BLAS = C3dB-new.F

  LIB_DEFINES = 


include ../../../config/makefile.h
ifeq ($(FC),pgf90)
 OBJ = $(OBJ_OPTIMIZE)
 OBJ_OPTIMIZE = 
endif
include ../../../config/makelib.h


ifeq ($(USE_MLIB),y)
  LIB_DEFINES += -DMLIB
endif

#ifeq ($(USE_MPIF),y)
#  LIB_DEFINES += -DMPI
#  LIB_INCLUDES += -I$(MPI_INCLUDE)
#else
  LIB_DEFINES += -DTCGMSG
#endif

