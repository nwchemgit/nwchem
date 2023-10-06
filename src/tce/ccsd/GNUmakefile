#$Id$

OBJ_OPTIMIZE = ccsd_e.o ccsd_t1.o ccsd_t2.o cc2_t1.o cc2_t2.o \
               ccsd_1prdm.o ccsd_1prdm_hh.o ccsd_1prdm_hp.o \
               ccsd_1prdm_ph.o ccsd_1prdm_pp.o \
               icsd_t1.o icsd_t2.o \
	       ccsd_kernels.o ccsd_t2_8.o tce_1b_dens_print.o

LIB_INCLUDES = -I../include

LIBRARY = libtce.a

USES_BLAS = ccsd_e.F ccsd_t1.F ccsd_t2.F cc2_t1.F cc2_t2.F \
            ccsd_1prdm_hh.F ccsd_1prdm_hp.F ccsd_1prdm_ph.F \
            ccsd_1prdm_pp.F ccsd_1prdm.F \
            icsd_t1.F icsd_t2.F ccsd_t2_8.F


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

