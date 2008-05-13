#$Id: GNUmakefile,v 1.1 2008-05-13 18:49:58 jhammond Exp $

OBJ_OPTIMIZE = \
ccsdt_lambda_o1.o \
ccsdt_lambda_o2.o \
ccsdt_lambda_o3.o \
ccsdt_y_tr1.o \
ccsdt_y_tr2.o \
ccsdt_y_tr3.o 

LIB_INCLUDES = -I../include

LIBRARY = libtce.a

USES_BLAS = \
ccsdt_lambda_o1.F \
ccsdt_lambda_o2.F \
ccsdt_lambda_o3.F \
ccsdt_y_tr1.F \
ccsdt_y_tr2.F \
ccsdt_y_tr3.F 
            
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

