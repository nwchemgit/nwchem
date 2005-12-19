#
# $Id: GNUmakefile,v 1.24 2005-12-19 21:24:32 marat Exp $
#


OBJ_OPTIMIZE =  qmmm.o task_qmmm_gradient.o mm_interface.o qmmm_link.o \
       qmmm_data.o qmmm_lbfgsb.o \
       qmmm_utils.o task_qmmm_energy.o qmmm_main.o \
       qmmm_check_forces.o qmmm_bq.o qmmm_esp.o qmmm_cons.o \
       qmmm_geom.o qmmm_bq_data.o \
       task_qmmm_optimize.o task_qmmm_dynamics.o qmmm_output.o \
       mm_utils.o qmmm_input.o

LIBRARY=libqmmm.a

HEADERS=

USES_BLAS=

include ../config/makefile.h

ifeq ($(TARGET),LINUX)
DEFINES += -DNEED_LOC
endif

include ../config/makelib.h

showobj:
	echo "obj $(OBJ) "
	echo "obj_opt $(OBJ_OPTIMIZE)"
