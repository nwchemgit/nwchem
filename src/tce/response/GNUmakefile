#$Id$

OBJ_OPTIMIZE = \
ccsd_o1.o \
ccsd_o2.o \
eomccsd_o1.o \
eomccsd_o2.o \
ccsdt_o1.o \
ccsdt_o2.o \
ccsdt_o3.o \
ccsdtq_o1.o \
ccsdtq_o2.o \
ccsdtq_o3.o \
ccsdtq_o4.o

LIB_INCLUDES = -I../include

LIBRARY = libtce.a

USES_BLAS = \
ccsd_o1.F \
ccsd_o2.F \
eomccsd_o1.F \
eomccsd_o2.F \
ccsdt_o1.F \
ccsdt_o2.F \
ccsdt_o3.F \
ccsdtq_o1.F \
ccsdtq_o2.F \
ccsdtq_o3.F \
ccsdtq_o4.F 
            
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

