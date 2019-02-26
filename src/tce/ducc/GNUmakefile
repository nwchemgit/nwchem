#$Id$

OBJ_OPTIMIZE =  ducc.o map_hcore.o rot_vir.o

LIB_INCLUDES = -I../include

LIBRARY = libtce.a

USES_BLAS = rot_vir.F
            
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

