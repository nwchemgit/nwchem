#
# $Id$
#

          OBJ = task_smd.o smd_input.o \
		chg_case.o ewald_self.o frce_shel.o \
		list_const.o md_inpt.o scle_velo.o \
		com_motion.o ewald_setp.o frce_shrt.o \
		list_excld.o  md_scle.o shll_qnch.o \
		erfc.o extrct_char.o init_velo.o \
		list_shell.o print_final.o tool_invrt.o \
		ewald_excl.o extrct_inte.o inte_leapf.o  \
		list_verlt.o print_output.o tool_rebox.o \
		ewald_real.o extrct_real.o inte_shake.o  \
		md_frce.o tool_randm.o tool_volme.o \
		ewald_recp.o frce_bond.o list_bonds.o  \
		md_init.o read_line.o verlt_test.o  \
                smd_init.o smd_end.o smd_param.o \
                smd_pdb.o  smd_atom.o smd_pot.o \
                smd_res.o smd_data.o smd_paramfile.o \
                task_smd_energy.o smd_rdinp.o smd_lat.o \
                smd_util.o smd_type.o smd_vel.o smd_force.o \
                smd_charge.o smd_coords.o smd_exlist.o \
                smd_system.o smd_ewald_data.o smd_ewald.o \
                smd_ewald_generic.o smd_cutoff_data.o \
                smd_vlist.o smd_pot_generic.o smd_energy.o \
                smd_shake.o smd_bond.o smd_shakelist.o \
                smd_bondlist.o smd_coulomb.o smd_rtdb_data.o \
                smd_coordfile.o task_smd_dynamics.o smd_system_data.o \
                smd_temper_data.o smd_mass.o smd_leapf.o \
                smd_leapf_shake.o smd_monomer.o smd_geom.o \
		process_task.o sched_groups.o smd_monomer_profile.o

      LIBRARY = libsmd.a


  LIB_DEFINES = -DMPI
ifneq ($(MPI_INCLUDE),)
  LIB_DEFINES += -I$(MPI_INCLUDE)
endif

include ../../config/makefile.h
include ../../config/makelib.h

