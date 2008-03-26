#$Id: GNUmakefile,v 1.2 2008-03-26 19:24:56 jhammond Exp $

OBJ_OPTIMIZE = beta_1_all.o beta_2_all.o beta_3_all.o beta_4_all.o beta_all_offsets.o \

LIB_INCLUDES = -I../include

LIBRARY = libtce.a

USES_BLAS = beta_1_all.F beta_2_all.F beta_3_all.F beta_4_all.F \
            
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

