#
# $Id: GNUmakefile,v 1.18 2005-05-27 17:59:21 marat Exp $
#

 OBJ =  mm_interface.o   
 OBJ_OPTIMIZE =  qmmm_utils.o task_qmmm_energy.o qmmm_main.o qmmm_data.o \
                 qmmm_geom.o  qmmm.o qmmm_input.o mm_utils.o qmmm_bq_data.o \
                 task_qmmm_gradient.o task_qmmm_optimize.o task_qmmm_dynamics.o qmmm_output.o \
                 qmmm_check_forces.o qmmm_bq.o qmmm_link.o 

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
