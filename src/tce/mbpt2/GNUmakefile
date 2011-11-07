#$Id$

OBJ_OPTIMIZE = mbpt2_e.o mbpt2_t1.o mbpt2_t2.o mbpt2_t2_b.o mbpt2_t3.o \
               mbpt2_t1_sd.o mbpt2_t2_sd.o

LIB_INCLUDES = -I../include

LIBRARY = libtce.a

USES_BLAS = mbpt2_e.F \
mbpt2_t1.F \
mbpt2_t2_b.F \
mbpt2_t2.F \
mbpt2_t3.F \
mbpt2_t1_sd.F \
mbpt2_t2_sd.F \
mbpt2_t4.F


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

