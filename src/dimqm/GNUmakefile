#
# $Id: GNUmakefile 20002 2011-02-26 02:05:02Z jochen $
#
#constants.o : constants.f90
#	      $(FC) -c constants.f90

   OBJ = dimqm_input.o dimpar_input.o dimqm_prep.o \
	 dimqm_main.o dimqm_polarizability.o dimqm_f2d.o dimqm_iterative.o \
	 qmr_real.o qmr_complex.o matvecReal.o matvecComplex.o \
         dimqm_EqmE.o t2.o \
         pphilo.o dimqm_fldpot.o dimqm_EqmN.o \
         dimqm_printAtomicDipoles.o dimqm_energy.o \
         dimqm_printDipole.o \
         dimqm_Dens2Dipoles.o dimqm_indDipoles.o \
         dimqm_printAtoms.o dimqm_addop.o dimqm_lclfld.o \
         seeded_qmr_real.o \
         qmr_seed_real.o seeded_qmr_real_augment.o \
         dimqm_tidy.o qmr_seed_complex.o seeded_qmr_complex.o \
         dimqm_broadcast.o CubicSpline.o dim_fock_xc.o \
         dim_getv.o dim_tabcd.o dim_eval_fnl.o \
         dim_grid_quadv0.o dim_grid_quadv0a.o dim_grid_quadv0b.o \
         dimqm_seeding.o dimqm_rohf.o dimqm_util.o dimqm_addop_uhf.o

   LIBRARY = libdimqm.a
#
# If you add -DJEFF_OUTPUT to LIB_DEFINES, the dipole polarizability will be printed in the same format as the TCE.
# This is very useful (1) if you wrote scripts to parse that particular formatting and (2) if you want the anisotrophy.
#
   LIB_DEFINES = -DDBLEPREC_RYS

   LIB_INCLUDES = -I../nwdft/include -I../ddscf -I../util -I../nwdft/grid -I../NWints/hondo

USES_BLAS = qmr_real.F dimqm_main.F 

HEADERS = dimqm.fh dimqm_constants.fh

include ../config/makefile.h
include ../config/makelib.h

#erf.o : erf.f90
#	$(FC) -c -fpp erf.f90
#FDEBUG += -O0
