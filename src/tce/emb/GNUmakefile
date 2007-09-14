#$Id: GNUmakefile,v 1.1 2007-09-14 19:04:01 kowalski Exp $

OBJ_OPTIMIZE = c1v.o c2vt1_s.o emb1.o emb2.o

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

