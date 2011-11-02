#$Id$

OBJ_OPTIMIZE = eomccsd_x1.o eomccsd_x2.o eomccsd_y1.o eomccsd_y2.o \
               eomccsd_density1.o eomccsd_denominator.o cc2_x1.o cc2_x2.o \
               eomccsd_1prdm.o eomccsd_1prdm_hh.o \
               eomccsd_1prdm_hp.o eomccsd_1prdm_ph.o eomccsd_1prdm_pp.o \
               eomccsd_x1_old.o eomccsd_x2_old.o \
               eomicsd_x1.o eomicsd_x2.o

LIB_INCLUDES = -I../include

LIBRARY = libtce.a

USES_BLAS = eomccsd_x1.F eomccsd_x2.F eomccsd_y1.F eomccsd_y2.F \
            eomccsd_denominator.F eomccsd_density1.F cc2_x1.F cc2_x2.F \
            eomccsd_1prdm_hh.F eomccsd_1prdm_hp.F \
            eomccsd_1prdm_ph.F eomccsd_1prdm_pp.F eomccsd_1prdm.F \
            eomccsd_x1_old.F eomccsd_x2_old.F \
            eomicsd_x1.F eomicsd_x2.F


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

