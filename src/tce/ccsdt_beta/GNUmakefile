#$Id: GNUmakefile,v 1.1 2008-05-02 04:28:28 jhammond Exp $

ifndef CCSDTLR
OBJ_OPTIMIZE = \
else
OBJ_OPTIMIZE = 

USES_BLAS = 

endif

LIB_INCLUDES = -I../include

LIBRARY = libtce.a

ifndef CCSDTLR
LIB_DEFINES = -DDEBUG_PRINT
else
LIB_DEFINES = -DDEBUG_PRINT -DCCSDTLR
endif
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

