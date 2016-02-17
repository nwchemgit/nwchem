#$Id$

OBJ_OPTIMIZE =  #tce_diis3.o tce_diis3c.o\
      tce_jacobi_ir1.o tce_jacobi_ir2.o tce_jacobi_ir3.o tce_jacobi_ir4.o \
      tce_jacobi2_ir1.o tce_jacobi2_ir2.o tce_jacobi2_ir3.o tce_jacobi2_ir4.o

LIB_INCLUDES = -I../include

LIBRARY = libtce.a

USES_BLAS = #tce_diis3.F tce_diis3c.F

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

