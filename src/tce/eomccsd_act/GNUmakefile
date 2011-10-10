#$Id$

OBJ_OPTIMIZE = eomccsd_x1a.o eomccsd_x2a.o q2na_hc1.o

LIB_INCLUDES = -I../include

LIBRARY = libtce.a

USES_BLAS = eomccsd_x1a.F eomccsd_x2a.F q2na_hc1.F


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

