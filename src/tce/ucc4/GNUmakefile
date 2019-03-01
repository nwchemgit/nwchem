#$Id$

OBJ_OPTIMIZE = ucc4_e_new.o  ucc4_t1_new.o  ucc4_t2.o  ucc4_t2_t2_v_t.o  ucc4_t3.o

LIB_INCLUDES = -I../include

LIBRARY = libtce.a

USES_BLAS = ucc4_e_new.F  ucc4_t1_new.F  ucc4_t2.F  ucc4_t2_t2_v_t.F  ucc4_t3.F

LIB_DEFINES = -DDEBUG_PRINT

# This replaces 3*TCE_SORT4+DGEMM with 6D loops (ccsd_kernels.F).
#LIB_DEFINES += -DUSE_LOOPS_NOT_DGEMM

#
# Possible #defines
#
#              -DDEBUG_PRINT          /* Debug printing routines */
#	       -DNOCOMMS              /* Turn off ALL communication calls */
#              -DBLOCK_TRANSF         /* Use block-transfer for comms calls */
#              -DOLD_AOBLOCK          /* Recover original AO blocking scheme */
#

HEADERS = 

include ../../config/makefile.h
include ../../config/makelib.h

