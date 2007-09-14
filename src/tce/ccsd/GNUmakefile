#$Id: GNUmakefile,v 1.3 2007-09-14 18:41:34 kowalski Exp $

OBJ_OPTIMIZE = ccsd_e.o ccsd_t1.o ccsd_t2.o cc2_t1.o cc2_t2.o

LIB_INCLUDES = -I../include

LIBRARY = libtce.a

USES_BLAS = ccsd_e.F ccsd_t1.F ccsd_t2.F cc2_t1.F cc2_t2.F

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

