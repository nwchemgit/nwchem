

OBJ = 

OBJ_OPTIMIZE =  mr_ccsd_e.o mr_ccsd_t1.o mr_ccsd_t2.o mr_ccsd_dcl.o \
		T2_C2_1.o T2_C2.o T2_pg_C2_1.o T2_pg_C2.o

LIB_INCLUDES = -I../../include -I../include

LIBRARY = libtce.a

USES_BLAS =	mr_ccsd_e.F mr_ccsd_t1.F mr_ccsd_t2.F mr_ccsd_dcl.F \
		T2_C2_1.F T2_C2.F T2_pg_C2_1.F T2_pg_C2.F

LIB_DEFINES = 


HEADERS = 

ifdef MRCC_METHODS
#ckbn mrcc-r-1 -3
      #LIB_DEFINES += -DMRCC_METHODS -DMRCC_LOCAL_FOCK -DMRCC_LOCAL_T1
      LIB_DEFINES += -DMRCC_METHODS #-DMRCC_LOCAL_FOCK -DMRCC_LOCAL_T1
endif

include ../../../config/makefile.h
include ../../../config/makelib.h

