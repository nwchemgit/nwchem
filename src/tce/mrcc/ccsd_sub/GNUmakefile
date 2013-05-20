
OBJ_OPTIMIZE =  ccsd_st1.o ccsd_st2.o ccsd_se.o ccsd_dcl_sub.o

LIB_INCLUDES = -I../../include -I../include

LIBRARY = libtce.a

USES_BLAS = ccsd_st1.F ccsd_st2.F ccsd_se.F ccsd_dcl_sub.F

#LIB_DEFINES += -DMRCC_METHODS -DMRCC_LOCAL_FOCK -DMRCC_LOCAL_T1
LIB_DEFINES += -DMRCC_METHODS #-DMRCC_LOCAL_FOCK -DMRCC_LOCAL_T1
#LIB_DEFINES = -DDEBUG_PRINT
#
# Possible #defines
#
#              -DDEBUG_PRINT          /* Debug printing routines */
#	       -DNOCOMMS              /* Turn off ALL communication calls */
#              -DBLOCK_TRANSF         /* Use block-transfer for comms calls */
#              -DOLD_AOBLOCK          /* Recover original AO blocking scheme */
#

HEADERS = 

include ../../../config/makefile.h
include ../../../config/makelib.h

