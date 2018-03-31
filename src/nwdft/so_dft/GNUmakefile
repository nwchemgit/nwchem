# $Id$
#
#      OBJ = dft_scf.o dft_s2_value.o dft_main0d.o
#      OBJ_OPTIMIZE = \
#	diis_driver.o	diis_bsolve.o	\
#	diis_core.o	diis_bld12.o diis_hamwgt.o diis_refvec.o	\
#	dft_dencvg.o dft_scfcvg.o	dft_densm.o	\
#	dft_invdiag.o	dft_init_guess.o \
#	dft_mxovl.o  dft_mxspin_ovlp.o
#

	OBJ = dft_scf_so.o \
	dft_densmso.o \
	dft_mxovlso.o \
	grad_force_so.o dft_gradients_so.o \
	vectors_so.o grad_dens_so.o grad1_so.o \
	diis_driver_so.o diis_bsolve_so.o diis_bld12_so.o\
	diis_core_so.o diis_hamwgt_so.o \
	diis_refvec_so.o  \
	movecs_anal_so.o \
	getv_coul.o getv_exch0.o getv_exchX.o getv_exchY.o getv_exchZ.o getv_dftxc.o\
	sym_mo_adapt_so.o sym_char_tab_so.o scf_sym_adap_so.o \
	sym_mo_ap_op_so.o sym_op_clsfy_so.o ga_orth_mos_so.o  \
	ga_orthog_so.o sym_op_cname_so.o grad_dens_so2.o \
	dft_scf_utils_so.o dft_frac_mo_so.o xc_exso.o getv_coul_cd.o
	OBJ_OPTIMIZE = 

     LIBRARY = libnwdft.a

 LIB_DEFINES =

LIB_INCLUDES = -I../include -I../../ddscf

   USES_BLAS = dft_mxovlso.F diis_bld12_so.F diis_bsolve_so.F \
sym_mo_adapt_so.F sym_op_clsfy_so.F dft_scf_so.F ga_orthog_so.F dft_scf_utils_so.F dft_densmso.F

include ../../config/makefile.h
include ../../config/makelib.h
