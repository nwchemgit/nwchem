#$Id$

OBJ_OPTIMIZE = mbpt3_e.o mbpt3_t1.o mbpt3_t2.o mbpt3_t1_sd.o mbpt3_t2_sd.o

LIB_INCLUDES = -I../include

LIBRARY = libtce.a

USES_BLAS = mbpt3_e.F \
mbpt3_t1.F \
mbpt3_t2.F \
mbpt3_t1_sd.F \
mbpt3_t2_sd.F


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

