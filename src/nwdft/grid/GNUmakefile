# $Id: GNUmakefile,v 1.16 2000-01-07 19:27:05 edo Exp $

      OBJ = \
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
        grid_qpts.o \
        grid_gen_pts.o \
        grid_file.o \
        grid_setspac.o \
        grid_beckew.o \
	grid_den.o \
	grid_den_2.o \
	grid_rchi.o \
	grid_rchi_2.o \
	grid_neighb.o \
	grid_signf.o \
	grid_ssw.o \
	grid_sym.o \
	grid_cleanup.o \
	grid_quadv0.o \
	grid_quadv0a.o \
        grid_balance.o \
        grid_gaussq.o \
        grid_roland.o \
        grid_bucket.o \
        grid_sphere.o \
	grid_utils.o

     LIBRARY = libnwdft.a

# LIB_DEFINES = -DDFT_TIMER
 LIB_DEFINES =

LIB_INCLUDES = -I../include

   USES_BLAS = grid_invariant.F grid_beckew.F grid_den_2.F grid_file.F grid_setspac.F grid_quadv0a.F



include ../../config/makefile.h
include ../../config/makelib.h


