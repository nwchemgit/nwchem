#$Id: GNUmakefile,v 1.1 2003-10-19 00:57:52 sohirata Exp $

OBJ_OPTIMIZE = eomccsdt_x1.o eomccsdt_x2.o eomccsdt_x3.o \
               eomccsdt_y1.o eomccsdt_y2.o eomccsdt_y3.o \
               eomccsdt_density1.o eomccsdt_denominator.o

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

