# $Id: GNUmakefile,v 1.80 2003-02-15 01:18:51 edo Exp $

	OBJ=xc_getv.o xc_exact_pot.o xc_sicutil.o xc_exso.o xc_chkgrad.o
	OBJ_OPTIMIZE = \
	a_dist.o \
	eumac.o \
	qdist.o \
	scat_mat.o \
	setACmat.o \
	xc_becke88.o \
	xc_beckew.o \
	xc_dirac.o \
	xc_dirac_screened.o \
	xc_fitv.o \
	xc_lebedev.o \
	xc_lyp88.o \
	xc_perdew86.o \
	xc_perdew91.o \
	xc_pw91lda.o \
	xc_eval_basis.o \
	xc_rhogen.o \
	xc_tabcd.o \
	xc_vwn_1.o \
	xc_vwn_1_rpa.o \
	xc_vwn_2.o \
	xc_vwn_3.o \
	xc_vwn_4.o \
	xc_vwn_5.o \
	xc_eval_fnl.o xc_fd.o xc_pbe96.o  \
	xc_hcth.o  xc_xpbe96.o xc_ssw.o  xc_xpw91.o xc_tarad.o xc_murakn.o \
	xc_signf.o  xc_gill96.o  xc_xmpw91.o\
	xc_oep.o xc_sic_pert.o xc_sic.o localize_sic.o \
	xc_optx.o xc_optc.o xc_atomblock.o xc_replicated.o xc_ft97.o


     LIBRARY = libnwdft.a

# LIB_DEFINES = -DDFT_TIMER
# LIB_DEFINES = -DQUADV0_REVERT

LIB_INCLUDES = -I../include


   USES_BLAS = xc_fitv.F xc_quadv0_a.F xc_quadv1.F xc_dirac.F xc_tabcd.F xc_rhogen.F xc_beckew.F xc_ssw.F xc_oep.F xc_sic_pert.F xc_sic.F localize_sic.F \
xc_eval_basis.F xc_fd.F  ../include/xc_rhogen.fh ../include/xc_tabcd.fh xc_optc.F xc_atomblock.F

include ../../config/makefile.h
include ../../config/makelib.h


