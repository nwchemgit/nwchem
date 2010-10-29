#$Id$

OBJ_OPTIMIZE = ccd_e.o ccd_t2.o

LIB_INCLUDES = -I../include

LIBRARY = libtce.a

USES_BLAS = ccd_e.F ccd_t2.F

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

