#ckbn only used if export EACCSD=y while building
 
OBJ_OPTIMIZE = eaccsd_x1.o eaccsd_x2.o tce_print_eax1.o tce_print_eax2.o tce_eax1_offset.o tce_eax2_offset.o \
               tce_eom_eaxguess.o tce_jacobi_eax1.o tce_jacobi_eax2.o

LIB_INCLUDES = -I../include

LIBRARY = libtce.a

USES_BLAS = eaccsd_x1.F eaccsd_x2.F tce_print_eax1.F tce_print_eax2.F tce_eax1_offset.F tce_eax2_offset.F \
            tce_eom_eaxguess.F tce_jacobi_eax1.F tce_jacobi_eax2.F


ifdef EACCSD
      LIB_DEFINES += -DEACCSD
endif
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

