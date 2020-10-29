#
# $Id$
#
      SUBDIRS = paw_utilities

 OBJ_OPTIMIZE = \
		nwpw_filefind.o \
		nwpw_scratch.o \
		nwpw_list.o \
		nwpw_scf_mixing.o \
		nwpw_spline.o \
		nwpw_interp.o \
		kerker_G.o \
		cpmd_properties.o \
		nwpw_rdf.o \
		pair_dist.o \
		pair_dist2.o \
		pair_distab.o \
		auto_corr.o \
		generate_unfolded_xyz.o \
		matrix.o simp.o \
		is_sORd.o \
		psi_filefind.o \
		v_psi_filefind.o \
		getfilling.o \
		CIF_write.o \
		densityofstates.o \
		MOTION_nframes.o \
		spdf_name.o \
		cell.o \
		cellgeometry.o\
		cpsi_data.o \
		psi_data.o \
		butter.o \
		print_elapsed_time.o \
		log_integrate.o \
		nwpw_parse.o \
		nwpw_matrix_invert.o \
		analysis_expansion_coef.o \
		nwpw_double_factorial.o \
		nwpw_kbpp_ray.o \
		nwpw_talker.o \
		nwpw_expression.o \
		nwpw_emachine.o \
		BA_memory.o btdb.o beom.o \
		task_paulsmatrix.o \
		nwpw_cosmo.o \
		nwpw_born.o \
		nwpw_cdft.o \
		nwpw_fftpack3d.o 

	  OBJ = 

      LIBRARY = libnwpwlib.a

    USES_BLAS = auto_corr.F  matrix.f nwpw_scf_mixing.F kerker_G.F cell.F cellgeometry.F nwpw_list.F cpsi_data.F psi_data.F nwpw_diis.F nwpw_kain.F nwpw_matrix_invert.F generate_unfolded_xyz.F nwpw_fftpack3d.F 

  LIB_DEFINES = 

   HEADERS = bafdecls.fh btdb.fh beom.fh



include ../../../config/makefile.h
include ../../../config/makelib.h

