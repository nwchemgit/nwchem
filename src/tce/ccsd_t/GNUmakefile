#$Id: GNUmakefile,v 1.1 2005-12-16 21:57:03 kowalski Exp $

OBJ_OPTIMIZE = ccsd_t.o ccsd_t_doubles.o ccsd_t_singles.o cr_ccsd_t.o cr_ccsd_t_D.o cr_ccsd_t_E.o cr_ccsd_t_N.o \
               ccsd2_t.o ccsd2_t_left.o ccsd2_q.o ccsd2_q_left.o ccsd2_q_right.o ccsdt2_q.o ccsdt2_q_left.o ccsdt2_q_right.o \
               lr_ccsd_t.o lr_ccsd_tq.o

LIB_INCLUDES = -I../include

LIBRARY = libtce.a

USES_BLAS = 

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

