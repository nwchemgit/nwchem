#$Id: GNUmakefile,v 1.2 2003-10-19 00:57:11 sohirata Exp $

OBJ_OPTIMIZE = eomccsd_x1.o eomccsd_x2.o eomccsd_y1.o eomccsd_y2.o eomccsd_density1.o eomccsd_denominator.o

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

