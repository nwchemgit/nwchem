# $Id: GNUmakefile,v 1.10 1994-09-02 01:08:58 gg502 Exp $

LIBRARY =	librimp2.a

OBJ =		rimp2_v_e2.o \
		mk_1d_shell_dist.o	ga_create_one_by_shells.o \
		xf3ci_step12.o	xf3ci_step3.o xform_3cints.o \
		occ_ranges.o	rimp2_xf_driver.o \
		close_3ci_file.o	open_3ci_file.o	put_3cints.o \
		get_3cints.o	banner.o	rimp2.o	\
		dgezer.o	dgecop.o	lnblnk.o \
		dgewr.o		atoi.o		indint.o \
		int_2c_ga.o     max_v_diagonal.o	rimp2_get_param.o \
		rimp2_timers.o  mk_fit_xf.o     drop_orbitals.o \
                mk_correl_map.o print_correl_map.o

OBJ_OPTIMIZE =  ga_matpow.o     rimp2_denom.o   rimp2_pair_energy.o \
                sum_pairs.o     xf3ci_step12a.o

USES_BLAS =	rimp2_v_e2.F    xf3ci_step12.F  xf3ci_step12a.F

include ../config/makefile.h
include ../config/makelib.h
