# $Id$

      OBJ_OPTIMIZE= \
	grid_pscreen.o \
	grid_sphere.o \
	grid_ssw.o \
	grid_sswold.o \
	grid_beckew.o \
	grid_signf.o \
	grid_quadvw.o \
	grid_quadvw1.o \
	grid_bucket.o \
	grid_utils.o  \
	grid_sympts.o  \
	grid_qpts.o  \
	grid_quadv0.o \
	grid_quadv0a.o \
	grid_quadv0b.o \
	grid_quadvc.o \
	grid_file.o grid_eaf.o
      OBJ= \
	grid_sic.o \
	grid_checks.o \
	grid_balance.o \
	grid_rtrunc.o \
	grid_atom_type_info.o \
	grid_setquad.o \
	grid_setrad.o \
	grid_setang.o \
	grid_invariant.o \
	grid_list.o \
	grid_setspac_params.o \
	grid_acc_def.o \
	grid_memory.o \
	grid_update_rtdb.o \
	grid_tasks.o \
	grid_gen_pts.o \
	grid_setspac.o \
	grid_neighb.o \
	grid_cleanup.o \
	grid_gaussq.o \
	grid_finpts.o \
	grid_init.o \
	grid_loop.o \
	grid_xcfit.o \
	dftg_loop.o

     LIBRARY = libnwdft.a

 LIB_DEFINES =

LIB_INCLUDES = -I../include

   USES_BLAS = grid_invariant.F grid_beckew.F grid_file.F grid_setspac.F \
 grid_quadv0a.F grid_quadv0b.F grid_quadvw.F grid_tasks.F grid_pscreen.F \
grid_quadvw1.F grid_xcfit.F grid_sympts.F grid_ssw.F grid_finpts.F grid_ylm.F grid_eaf.F

include ../../config/makefile.h
include ../../config/makelib.h
