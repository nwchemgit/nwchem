#$Id: GNUmakefile,v 1.3 2005-12-20 21:19:40 edo Exp $

OBJ_OPTIMIZE = cisdtq_e.o cisdtq_c1.o cisdtq_c2.o cisdtq_c3.o cisdtq_c4.o

LIB_INCLUDES = -I../include

LIBRARY = libtce.a

USES_BLAS = cisdtq_c1.F \
cisdtq_c2.F \
cisdtq_c3.F \
cisdtq_c4.F \
cisdtq_e.F


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

