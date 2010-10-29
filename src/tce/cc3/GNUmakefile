#$Id$

OBJ_OPTIMIZE = cc3_t3.o ccsdt1b_t3.o ccsdt2_t3.o ccsdt3_t3.o eomccsdt1_x3.o eomcc3_x3.o eomccsdt3_x3.o cc3_lambda3.o

LIB_INCLUDES = -I../include

LIBRARY = libtce.a

USES_BLAS = cc3_t3.F ccsdt1b_t3.F ccsdt2_t3.F ccsdt3_t3.F eomccsdt1_x3.F eomcc3_x3.F eomccsdt3_x3.F cc3_lambda3.F

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

