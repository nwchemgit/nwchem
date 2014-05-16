#ckbn only used if export IPCCSD=y while building

OBJ_OPTIMIZE = ipccsd_x1.o ipccsd_x2.o tce_print_ipx1.o tce_print_ipx2.o tce_ipx1_offset.o tce_ipx2_offset.o \
               tce_eom_ipxguess.o tce_jacobi_ipx1.o tce_jacobi_ipx2.o 

LIB_INCLUDES = -I../include

LIBRARY = libtce.a

USES_BLAS = ipccsd_x1.F ipccsd_x2.F tce_print_ipx1.F tce_print_ipx2.F tce_ipx1_offset.F tce_ipx2_offset.F \
            tce_eom_ipxguess.F tce_jacobi_ipx1.F tce_jacobi_ipx2.F 


ifdef IPCCSD
      LIB_DEFINES += -DIPCCSD
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

