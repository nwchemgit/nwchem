#$Id: GNUmakefile,v 1.1 2003-09-16 21:54:12 sohirata Exp $

OBJ_OPTIMIZE = lccsd_e.o lccsd_t1.o lccsd_t2.o

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

