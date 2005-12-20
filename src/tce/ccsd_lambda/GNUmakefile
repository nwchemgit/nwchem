#$Id: GNUmakefile,v 1.2 2005-12-20 21:19:40 edo Exp $

OBJ_OPTIMIZE = ccsd_lambda1.o ccsd_lambda2.o ccsd_density1.o

LIB_INCLUDES = -I../include

LIBRARY = libtce.a

USES_BLAS = ccsd_density1.F \
ccsd_lambda1.F \
ccsd_lambda2.F


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

