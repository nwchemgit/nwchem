#
# $Id: GNUmakefile,v 1.15 2004-12-21 02:15:21 marat Exp $
#

 OBJ_OPTIMIZE =  qmmm_utils.o task_qmmm_energy.o qmmm_main.o qmmm_data.o qmmm_bq_data.o mm_interface.o  \
                 qmmm_geom.o  qmmm.o qmmm_input.o mm_utils.o \
                 task_qmmm_gradient.o task_qmmm_optimize.o task_qmmm_dynamics.o qmmm_output.o \
                 qmmm_check_forces.o

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
