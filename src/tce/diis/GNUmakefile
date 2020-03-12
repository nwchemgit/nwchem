#$Id$

OBJ_OPTIMIZE = tce_diis.o tce_diis_init.o tce_diis_switch.o tce_diis_tidy.o \
      tce_jacobi_t1.o tce_jacobi_t2.o tce_jacobi_t3.o tce_jacobi_t4.o \
      tce_jacobi_l1.o tce_jacobi_l2.o tce_jacobi_l3.o tce_jacobi_l4.o \
      tce_jacobi_x1.o tce_jacobi_x2.o tce_jacobi_x3.o tce_jacobi_x4.o \
      tce_jacobi_y1.o tce_jacobi_y2.o tce_jacobi_y3.o tce_jacobi_y4.o \
      tce_jacobi_t3a.o tce_jacobi_x3a.o t2_2_y2.o

LIB_INCLUDES = -I../include

LIBRARY = libtce.a

USES_BLAS = tce_diis.F
            
LIB_DEFINES = -DDEBUG_PRINT

CCSDT=1
ifdef CCSDT
    LIB_DEFINES += -DCCSDT
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

