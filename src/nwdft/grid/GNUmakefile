# $Id: GNUmakefile,v 1.6 1999-06-14 22:39:27 d3h449 Exp $

      OBJ = \
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
	grid_balance.o \
	grid_utils.o

     LIBRARY = libnwdft.a

# LIB_DEFINES = -DDFT_TIMER
 LIB_DEFINES =

LIB_INCLUDES = -I../include

   USES_BLAS = grid_invariant.F

include ../../config/makefile.h
include ../../config/makelib.h


