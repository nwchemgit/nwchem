#$Id$

OBJ_OPTIMIZE = cisdt_e.o cisdt_c1.o cisdt_c2.o cisdt_c3.o

LIB_INCLUDES = -I../include

LIBRARY = libtce.a

USES_BLAS = cisdt_c1.F \
cisdt_c2.F \
cisdt_c3.F \
cisdt_e.F


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

