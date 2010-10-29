#$Id$

OBJ_OPTIMIZE = ccsd_lambda1_1.o ccsd_lambda1_2.o ccsd_lambda2.o ccsd_density1.o ccsd_lambda1_mod.o

LIB_INCLUDES = -I../include

LIBRARY = libtce.a

USES_BLAS = ccsd_density1.F \
ccsd_lambda1_1.F ccsd_lambda1_2.F \
ccsd_lambda2.F ccsd_lambda1_mod.F


LIB_DEFINES = -DDEBUG_PRINT
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

