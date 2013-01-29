#
# $Id$
#
OBJ_OPTIMIZE =  qmd_init.o qmd_driver.o qmd_main.o qmd_input.o\
		qmd_start.o task_qmd.o qmd_nucl_dt2.o qmd_nucl_dt.o\
		qmd_gradient.o qmd_geom_store.o qmd_thermostat.o\
		qmd_kinetic.o qmd_write_trj.o qmd_write_rst.o

LIBRARY=libqmd.a

HEADERS=

USES_BLAS=

include ../config/makefile.h
include ../config/makelib.h

showobj:
	echo "obj $(OBJ) "
	echo "obj_opt $(OBJ_OPTIMIZE)"
