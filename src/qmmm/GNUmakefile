#
# $Id: GNUmakefile,v 1.14 2004-12-17 02:21:48 marat Exp $
#

 OBJ_OPTIMIZE =  qmmm_utils.o task_qmmm_energy.o qmmm_main.o qmmm_data.o qmmm_bq_data.o mm_interface.o  \
                 qmmm_geom.o  qmmm.o qmmm_input.o mm_utils.o \
                 task_qmmm_gradient.o task_qmmm_optimize.o task_qmmm_dynamics.o qmmm_output.o

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
