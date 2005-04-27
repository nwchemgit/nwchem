#
# $Id: GNUmakefile,v 1.17 2005-04-27 00:55:46 marat Exp $
#

 OBJ =  mm_interface.o   
 OBJ_OPTIMIZE =  qmmm_utils.o task_qmmm_energy.o qmmm_main.o qmmm_data.o \
                 qmmm_geom.o  qmmm.o qmmm_input.o mm_utils.o qmmm_bq_data.o \
                 task_qmmm_gradient.o task_qmmm_optimize.o task_qmmm_dynamics.o qmmm_output.o \
                 qmmm_check_forces.o qmmm_bq.o

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
