#$Id: GNUmakefile,v 1.1 2007-08-11 01:21:37 jhammond Exp $

OBJ_OPTIMIZE = \
ccsd_o1.o \
second_ccsd_o2.o \
ccsdt_o1.o \
second_ccsdt_o2.o \
second_ccsdt_o3.o \
ccsdtq_o1.o \
second_ccsdtq_o2.o \
second_ccsdtq_o3.o \
third_ccsdtq_o4.o

LIB_INCLUDES = -I../include

LIBRARY = libtce.a

USES_BLAS = \
ccsd_o1.F \
second_ccsd_o2.F \
ccsdt_o1.F \
second_ccsdt_o2.F \
second_ccsdt_o3.F \
ccsdtq_o1.F \
second_ccsdtq_o2.F \
second_ccsdtq_o3.F \
third_ccsdtq_o4.F 
            
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

