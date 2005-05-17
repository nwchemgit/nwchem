#
# $Id: GNUmakefile,v 1.24 2005-05-17 01:55:56 bylaska Exp $
#

 OBJ_OPTIMIZE = \
		nwpw_filefind.o \
		nwpw_scratch.o \
		nwpw_list.o \
		nwpw_scf_mixing.o \
		kerker_G.o \
		cpmd_properties.o \
		pair_dist.o \
		auto_corr.o \
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
		cellgeometry.o

	  OBJ = 

      LIBRARY = libnwpwlib.a

    USES_BLAS = auto_corr.F  matrix.f nwpw_scf_mixing.F kerker_G.F cell.F

  LIB_DEFINES = 


include ../../../config/makefile.h
include ../../../config/makelib.h

