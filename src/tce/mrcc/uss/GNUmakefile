

OBJ = 

OBJ_OPTIMIZE =  tce_uss_offdiagonal_1.o tce_uss_offdiagonal_2.o \
               tce_uss_offdiagonal_3.o tce_uss_offdiagonal_4.o \
               ccsd_r3.o ccsd_r4.o tce_uss_aux.o tce_uss_T1_R1.o tce_uss_T2_R2.o

LIB_INCLUDES = -I../../include -I../include

LIBRARY = libtce.a

USES_BLAS =	tce_uss_offdiagonal_1.F tce_uss_offdiagonal_2.F \
                tce_uss_offdiagonal_3.F tce_uss_offdiagonal_4.F \
                ccsd_r3.F ccsd_r4.F tce_uss_aux.F tce_uss_T1_R1.F tce_uss_T2_R2.F

LIB_DEFINES = 


HEADERS = 

ifdef MRCC_METHODS
#ckbn mrcc-r-1 -3
      #LIB_DEFINES += -DMRCC_METHODS -DMRCC_LOCAL_FOCK -DMRCC_LOCAL_T1
      LIB_DEFINES += -DMRCC_METHODS #-DMRCC_LOCAL_FOCK -DMRCC_LOCAL_T1
endif

include ../../../config/makefile.h
include ../../../config/makelib.h

