#
# $Id: GNUmakefile,v 1.12 2008-05-14 18:14:29 d3p708 Exp $
#

 OBJ_OPTIMIZE = C3dB-new.o c_trans.o C3dB_pfft.o
	  OBJ = 

ifeq ($(USE_MPIF),y)
  OBJ_OPTIMIZE += C3dB-mpi.o
else
  OBJ_OPTIMIZE += C3dB-tcgmsg.o
endif



      LIBRARY = libnwpwlib.a

    USES_BLAS = C3dB-new.F C3dB_pfft.F

  LIB_DEFINES = 


include ../../../config/makefile.h
#ifeq ($(FC),pgf90)
 #OBJ = C3dB-new.o
 #OBJ_OPTIMIZE = c_trans.o
#endif
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


