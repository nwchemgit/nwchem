

OBJ = 

OBJ_OPTIMIZE = mr_ccsd_e.o mr_ccsd_t1.o mr_ccsd_t2.o mr_ccsd_dcl.o

LIB_INCLUDES = -I../../include -I../include

LIBRARY = libtce.a

USES_BLAS = 

LIB_DEFINES = 


HEADERS = 

ifdef MRCC_METHODS
#ckbn mrcc-r-1 -3
      LIB_DEFINES += -DMRCC_METHODS
endif

include ../../../config/makefile.h
include ../../../config/makelib.h

