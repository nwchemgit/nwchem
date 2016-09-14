#
# $Id$
#
OBJ_OPTIMIZE =  qmd_init.o qmd_driver.o qmd_main.o qmd_input.o\
		qmd_start.o task_qmd.o qmd_nucl_dt2.o qmd_nucl_dt.o\
		qmd_gradient.o qmd_geom_store.o qmd_thermostat.o\
		qmd_kinetic.o qmd_write_trj.o qmd_write_rst.o\
		qmd_ran1.o qmd_gasdev.o qmd_randomize.o\
                qmd_nucl_pos.o qmd_nucl_vel.o qmd_read_rst.o\
                qmd_gamdev.o qmd_sumno.o qmd_fix_momentum.o\
                qmd_freeze_inactive.o qmd_cavity.o\
		qmd_energy.o qmd_coup.o qmd_move_civecs.o\
		qmd_elec_prop.o qmd_tdks_enmo.o qmd_tddft_enmo.o\
		qmd_namd_init.o qmd_namd_input.o qmd_set_den.o\
		qmd_coords_vels.o qmd_civecs_read.o qmd_property.o

LIBRARY=libqmd.a

HEADERS= ../util/util_params.fh

USES_BLAS=

include ../config/makefile.h
include ../config/makelib.h

showobj:
	echo "obj $(OBJ) "
	echo "obj_opt $(OBJ_OPTIMIZE)"
