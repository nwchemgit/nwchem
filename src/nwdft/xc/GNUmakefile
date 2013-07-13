# $Id$

HEADERS = xc.fh xc_vdw.fh xc_params.fh

	OBJ = xc_getv.o \
	xc_exact_pot.o \
	xc_sicutil.o \
	xc_exso.o \
	chk_dens.o \
	xc_util.o \
	xc_chkdispauto.o \
	xc_vdw_init.o 2ndDerivB97.o 

	OBJ_OPTIMIZE = xc_eval_basis.o  \
	a_dist.o \
	eumac.o \
	qdist.o \
	scat_mat.o \
	setACmat.o \
	xc_becke88.o \
	xc_dirac.o \
	xc_bnl.o \
	xc_dirac_screened.o \
	xc_fitv.o \
	xc_hfexch.o\
	xc_lebedev.o \
	xc_lyp88.o \
	xc_perdew86.o \
	xc_perdew91.o \
	xc_pw91lda.o xc_pw91ldag.o \
	xc_rhogen.o \
	xc_tabcd.o \
	xc_vwn_1.o \
	xc_vwn_1_rpa.o \
	xc_vwn_2.o \
	xc_vwn_3.o \
	xc_vwn_4.o \
	xc_vwn_5.o \
	xc_eval_fnl.o xc_fd.o xc_pbe96.o  \
	xc_hcth.o  xc_xpbe96.o xc_xpw91.o xc_tarad.o xc_murakn.o \
	xc_signf.o  xc_gill96.o  xc_xmpw91.o  xc_xpkzb99.o \
	xc_cpkzb99.o   xc_cMpbe96.o xc_xtpss03.o xc_ctpss03.o \
	xc_oep.o xc_sic_pert.o xc_sic.o localize_sic.o \
	xc_optx.o xc_optc.o xc_atomblock.o xc_replicated.o xc_ft97.o \
	xc_lb94.o xc_cs00.o xc_op.o xc_kop.o \
	xc_bc95.o xc_xm05.o xc_cm05.o xc_xpw6.o\
	xc_xvs98.o xc_cvs98.o xc_xm06.o xc_cm06.o \
	xc_spbe96.o xc_kt1.o xc_ssbD_1.o xc_s12x.o xc_cams12x.o \
	xc_camb88.o xc_att_xc.o xc_camxpbe96.o xc_camxlsd.o xc_xwpbe.o hse08fx.o \
	xc_vdw_main.o xc_vdw_util.o \
	xc_xsogga.o \
	xc_xm11.o \
	xc_csogga.o\
	xc_cm11.o \
	steric_energy.o

     LIBRARY = libnwdft.a

# LIB_DEFINES = -DQUADV0_REVERT

LIB_INCLUDES = -I../include
ifeq ($(TARGET),cray-sv2)
    EXPLICITF=TRUE
endif

   USES_BLAS = xc_fitv.F xc_dirac.F xc_bnl.F xc_att_xc.F xc_tabcd.F xc_rhogen.F \
               xc_oep.F xc_sic_pert.F xc_sic.F localize_sic.F \
               xc_eval_basis.F xc_fd.F xc_optc.F xc_atomblock.F scat_mat.F \
               xc_replicated.F xc_vdw_main.F xc_vdw_util.F

include ../../config/makefile.h
include ../../config/makelib.h
