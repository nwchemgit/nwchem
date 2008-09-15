
ifdef CCSDTLR
OBJ_OPTIMIZE = ccsdt_lr_alpha1.o ccsdt_lr_alpha2.o
USES_BLAS =  ccsdt_lr_alpha1.F ccsdt_lr_alpha2.F
endif

LIB_INCLUDES = -I../include

LIBRARY = libtce.a

ifndef CCSDTLR
LIB_DEFINES = -DDEBUG_PRINT
else
LIB_DEFINES = -DDEBUG_PRINT -DCCSDTLR
endif
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

