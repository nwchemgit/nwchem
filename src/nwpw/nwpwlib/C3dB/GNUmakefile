#
# $Id$
#

 OBJ_OPTIMIZE = C3dB-new.o c_trans.o C3dB_pfft.o
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
  OBJ_OPTIMIZE += C3dB-mpi.o
else
  OBJ_OPTIMIZE += C3dB-tcgmsg.o
endif



      LIBRARY = libnwpwlib.a

    USES_BLAS = C3dB-new.F C3dB_pfft.F  C3dB-tcgmsg.F C3dB-mpi.F

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
else
  LIB_DEFINES += -DTCGMSG
endif

ifeq ($(USE_MPIF4),y)
  LIB_DEFINES += -DMPI4
endif


