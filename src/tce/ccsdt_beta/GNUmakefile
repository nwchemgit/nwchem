#$Id: GNUmakefile,v 1.2 2008-05-05 02:26:58 jhammond Exp $

#ifndef CCSDTLR
#OBJ_OPTIMIZE = ccsdt_lr_beta_1.o ccsdt_lr_beta_2.o ccsdt_lr_beta_3.o ccsdt_lr_beta_4.o
#else
OBJ_OPTIMIZE = ccsdt_lr_beta_1.o ccsdt_lr_beta_2.o ccsdt_lr_beta_3.o ccsdt_lr_beta_4.o

USES_BLAS = ccsdt_lr_beta_1.F ccsdt_lr_beta_2.F ccsdt_lr_beta_3.F ccsdt_lr_beta_4.F

#endif

LIB_INCLUDES = -I../include

LIBRARY = libtce.a

#ifndef CCSDTLR
#LIB_DEFINES = -DDEBUG_PRINT
#else
LIB_DEFINES = -DDEBUG_PRINT -DCCSDTLR
#endif
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

