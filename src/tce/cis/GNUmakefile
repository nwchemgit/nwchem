#$Id: GNUmakefile,v 1.2 2005/12/20 21:19:40 edo Exp $

OBJ_OPTIMIZE =  cis.o cis_x1.o cis_xguess.o cis_xdiagon.o cis_xtidy.o

LIB_INCLUDES = -I../include

LIBRARY = libtce.a

USES_BLAS = cis_x1.F cis_xdiagon.F \


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

