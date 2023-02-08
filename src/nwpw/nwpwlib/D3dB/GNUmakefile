#
# $Id$
#

 OBJ_OPTIMIZE = D3dB-new.o D3dBs-new.o generate_map_indexes.o hilbert.o hcurve.o olist.o trans.o D3dB_pfft.o D3dBs_pfft.o
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
  OBJ_OPTIMIZE += D3dB-mpi.o D3dBs-mpi.o
else
  OBJ_OPTIMIZE += D3dB-tcgmsg.o
endif


      LIBRARY = libnwpwlib.a

    USES_BLAS = D3dB-new.F D3dB_pfft.F D3dB-mpi.F D3dB-tcgmsg.F D3dBs-new.F D3dBs_pfft.F D3dBs-mpi.F
  LIB_DEFINES = 


include ../../../config/makefile.h
include ../../../config/makelib.h


ifeq ($(USE_MLIB),y)
  LIB_DEFINES += -DMLIB
endif

ifeq ($(USE_MPIF),y)
  LIB_DEFINES += -DMPI
else
  USES_BLAS += D3dB-tcgmsg.F
  LIB_DEFINES += -DTCGMSG
endif

ifeq ($(USE_MPIF4),y)
LIB_DEFINES += -DMPI4
endif

ifeq ($(USE_FFTW3),y)
  LIB_DEFINES += -DFFTW3
  LIB_INCLUDES += -I$(FFTW3_INCLUDE)
endif
