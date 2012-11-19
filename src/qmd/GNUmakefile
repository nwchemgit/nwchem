#
# $Id$
#
OBJ_OPTIMIZE =  qmd_init.o qmd_driver.o qmd_main.o qmd_input.o \
		qmd_start.o task_qmd.o qmd_nucl_step.o \
		qmd_veloverlet.o qmd_kinetic.o qmd_berendsen.o qmd_nosehoover.o

LIBRARY=libqmd.a

HEADERS=

USES_BLAS=

include ../config/makefile.h
include ../config/makelib.h

showobj:
	echo "obj $(OBJ) "
	echo "obj_opt $(OBJ_OPTIMIZE)"
