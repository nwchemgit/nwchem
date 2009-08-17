#$Id: GNUmakefile,v 1.11 2008-10-16 14:30:28 jhammond Exp $

OBJ          = cr_ccsd_t.o lr_ccsd_t.o

OBJ_OPTIMIZE = ccsd_t.o ccsd_t_doubles.o ccsd_t_singles.o lambda_ccsd_t.o lambda_ccsd_t_left.o\
               ccsd_t_doubles_l.o ccsd_t_singles_l.o \
               ccsd2_t.o ccsd2_t_left.o ccsd2_q.o ccsd2_q_left.o ccsd2_q_right.o \
               ccsdt2_q.o ccsdt2_q_left.o ccsdt2_q_right.o \
               cr_ccsd_t_D.o cr_ccsd_t_E.o cr_ccsd_t_N.o lr_ccsd_tq.o \
               cr_qcisd_t.o m2_dot_m1t1.o ccsd_t_restart.o ccsd_t_sse2.o

LIB_INCLUDES = -I../include

LIBRARY = libtce.a

USES_BLAS = ccsd_t_doubles.F ccsd_t_singles.F \
            ccsd_t_doubles_l.F ccsd_t_singles_l.F \
            lambda_ccsd_t_left.F \
            ccsd2_q_left.F ccsd2_q_right.F ccsd2_t_left.F \
            ccsdt2_q_left.F ccsdt2_q_right.F \
            cr_ccsd_t_D.F cr_ccsd_t_E.F cr_ccsd_t_N.F lr_ccsd_tq.F \
            cr_qcisd_t.F m2_dot_m1t1.F


LIB_DEFINES += -DDEBUG_PRINT
#LIB_DEFINES += -DSSE2_VERSION
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

