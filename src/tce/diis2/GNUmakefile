#$Id$

OBJ_OPTIMIZE = tce_diis2.o \
      tce_jacobi_tr1.o tce_jacobi_tr2.o tce_jacobi_tr3.o tce_jacobi_tr4.o \
      tce_jacobi_lr1.o tce_jacobi_lr2.o tce_jacobi_lr3.o tce_jacobi_lr4.o

LIB_INCLUDES = -I../include

LIBRARY = libtce.a

USES_BLAS = tce_diis2.F
            
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

