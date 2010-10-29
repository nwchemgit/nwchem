#$Id$

OBJ_OPTIMIZE = \
ccsd_y_o1.o \
ccsd_y_o2.o \
ccsd_yr1.o \
ccsd_yr2.o \
ccsd_y_tr1.o \
ccsd_y_tr2.o \
alpha_2_yr.o \
ccsd_lr_d1_hh.o \
ccsd_lr_d1_ph.o \
ccsd_lr_d1_hp.o \
ccsd_lr_d1_pp.o

LIB_INCLUDES = -I../include

LIBRARY = libtce.a

USES_BLAS = \
ccsd_y_o1.F \
ccsd_y_o2.F \
ccsd_yr1.F \
ccsd_yr2.F \
ccsd_y_tr1.F \
ccsd_y_tr2.F \
alpha_2_yr.F \
ccsd_lr_d1_hh.F \
ccsd_lr_d1_ph.F \
ccsd_lr_d1_hp.F \
ccsd_lr_d1_pp.F
            
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

