#
# $Id: GNUmakefile,v 1.7 2004-04-20 17:52:52 bylaska Exp $
#

 OBJ_OPTIMIZE = D3dB-new.o generate_map_indexes.o hilbert.o olist.o trans.o
 OBJ = 

ifeq ($(USE_MPI),y)
  OBJ_OPTIMIZE += D3dB-mpi.o
else
  OBJ_OPTIMIZE += D3dB-tcgmsg.o
endif

      LIBRARY = libnwpwlib.a

    USES_BLAS = D3dB-new.F 

  LIB_DEFINES = 


include ../../../config/makefile.h
ifeq ($(FC),pgf90)
 OBJ_OPTIMIZE = 
 OBJ = D3dB-new.o
endif
include ../../../config/makelib.h


ifeq ($(USE_MLIB),y)
  FOPTIONS += -DMLIB
endif

ifeq ($(USE_MPI),y)
  FOPTIONS += -DMPI
  LIB_INCLUDES += -I$(MPI_INCLUDE)
else
  FOPTIONS += -DTCGMSG
endif

