#$Id: GNUmakefile,v 1.2 2005-12-16 22:45:45 kowalski Exp $

OBJ_OPTIMIZE = cisdtq_e.o cisdtq_c1.o cisdtq_c2.o cisdtq_c3.o cisdtq_c4.o

LIB_INCLUDES = -I../include

LIBRARY = libtce.a

USES_BLAS = 

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

