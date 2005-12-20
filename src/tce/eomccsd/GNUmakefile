#$Id: GNUmakefile,v 1.3 2005-12-20 21:19:40 edo Exp $

OBJ_OPTIMIZE = eomccsd_x1.o eomccsd_x2.o eomccsd_y1.o eomccsd_y2.o eomccsd_density1.o eomccsd_denominator.o

LIB_INCLUDES = -I../include

LIBRARY = libtce.a

USES_BLAS = eomccsd_denominator.F \
eomccsd_density1.F \
eomccsd_x1.F \
eomccsd_x2.F \
eomccsd_y1.F \
eomccsd_y2.F


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

