#
# $Id$
#


OBJ =  qmmm.o task_qmmm_gradient.o mm_interface.o qmmm_link.o \
       qmmm_data.o qmmm_lbfgsb.o \
       qmmm_utils.o task_qmmm_energy.o qmmm_main.o \
       qmmm_check_forces.o qmmm_bq.o qmmm_esp.o qmmm_cons.o \
       qmmm_geom.o qmmm_bq_data.o \
       task_qmmm_optimize.o task_qmmm_dynamics.o qmmm_output.o \
       mm_utils.o qmmm_input.o qmmm_pes.o qmmm_property.o task_qmmm_fep.o \
       task_qmmm_abc.o

LIBRARY=libqmmm.a

HEADERS=

USES_BLAS= task_qmmm_gradient.F qmmm_lbfgsb.F mm_interface.F 

include ../config/makefile.h

ifeq ($(TARGET),LINUX)
DEFINES += -DNEED_LOC
endif

include ../config/makelib.h

showobj:
	echo "obj $(OBJ) "
	echo "obj_opt $(OBJ_OPTIMIZE)"
