#$Id: GNUmakefile,v 1.1 2007-10-30 16:48:23 jhammond Exp $

OBJ_OPTIMIZE = ccsd_quad_ot1.o ccsd_quad_ot2.o \
               ccsd_quad_xy1.o ccsd_quad_xy2_1.o \
               ccsd_quad_xy2_2.o  ccsd_quad_xy2_3.o

LIB_INCLUDES = -I../include

LIBRARY = libtce.a

USES_BLAS = ccsd_quad_ot1.F ccsd_quad_ot2.F \
            ccsd_quad_xy1.F ccsd_quad_xy2_1.F \
            ccsd_quad_xy2_2.F  ccsd_quad_xy2_3.F

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

