#$Id: GNUmakefile,v 1.4 2007/09/19 22:18:50 d3p975 Exp $

OBJ_OPTIMIZE = ccsd_ea.o ccsd_t1a.o ccsd_t2a.o activity.o

LIB_INCLUDES = -I../include

LIBRARY = libtce.a

USES_BLAS = ccsd_ea.F ccsd_t1a.F ccsd_t2a.F

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

