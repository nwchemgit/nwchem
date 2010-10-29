#$Id$

OBJ_OPTIMIZE = cisd_e.o cisd_c1.o cisd_c2.o

LIB_INCLUDES = -I../include

LIBRARY = libtce.a

USES_BLAS = cisd_c1.F \
cisd_c2.F \
cisd_e.F


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

