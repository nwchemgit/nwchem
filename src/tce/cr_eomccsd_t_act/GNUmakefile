#$Id$

OBJ_OPTIMIZE = cr_ccsd_t_E_act.o cr_eomccsd_t_act.o q3rexpt2_act.o  \
               cr_ccsd_t_N_act.o creomccsd_t_n2_act_mem.o  \
               cr_eomccsd_tasim.o cr_eomt_simac.o ccsd_t_d_act.o

LIB_INCLUDES = -I../include

LIBRARY = libtce.a

USES_BLAS = cr_ccsd_t_E_act.F  q3rexpt2_act.F  \
            cr_ccsd_t_N_act.F creomccsd_t_n2_act_mem.F \
            cr_eomt_simac.F ccsd_t_d_act.F  cr_eomccsd_t_act.F \
            cr_eomccsd_tasim.F


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

