# $Id$

HEADERS = xc.fh xc_vdw.fh xc_params.fh xc_hcth.fh

	OBJ = xc_getv.o \
	xc_getvxc.o \
	xc_exact_pot.o \
	xc_sicutil.o \
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
	xc_hcth.o  xc_xpbe96.o xc_xpw91.o xc_tarad.o xc_murakn.o xc_mitani.o \
	xc_mitani_dmr.o \
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
	xc_xsogga.o xc_xmvs15.o \
	xc_xm11.o xc_xmn12.o xc_xmn15.o \
	xc_csogga.o\
	xc_cm11.o \
	steric_energy.o\
	xc_3rd_deriv.o xc_xpw86.o xc_xb86b.o xc_xdm.o \
	ts_eval_fnl.o ts_tf.o ts_vw.o xc_xn12.o \
	xc_xscan.o xc_cscan.o xc_xscanl.o xc_cscanl.o ts_pc.o\
	xc_xncap.o xc_xncap.o xc_xr2scan.o xc_xr2scanl.o\
	xc_cr2scan.o xc_cr2scanl.o xc_hcth_extra.o \
	xc_rhogen_adft.o xc_tabcd_adft.o

     LIBRARY = libnwdft.a

# LIB_DEFINES = -DQUADV0_REVERT


LIB_INCLUDES = -I../include

   USES_BLAS = xc_fitv.F xc_dirac.F xc_bnl.F xc_att_xc.F xc_tabcd.F xc_rhogen.F \
               xc_oep.F xc_sic_pert.F xc_sic.F localize_sic.F \
               xc_eval_basis.F xc_eval_fnl.F xc_fd.F xc_optc.F xc_atomblock.F scat_mat.F \
               xc_replicated.F xc_vdw_main.F xc_vdw_util.F\
               xc_3rd_deriv.F xc_hcth.F ts_tf.F xc_xn12.F xc_util.F \
	       xc_tabcd_adft.F xc_rhogen_adft.F

include ../../config/makefile.h
include ../../config/makelib.h
ifdef EXTERNAL_GA_PATH
  GA_WAPI.H=${EXTERNAL_GA_PATH}/include/ga-wapi.h
else
  GA_WAPI.H=${NWCHEM_TOP}/src/tools/build/global/src/ga-wapi.h
endif
HAVE_SET_GA_PROPERTY = $(shell cat ${GA_WAPI.H} | awk ' /wnga_set_property/ {print "Y"}')
ifeq ($(HAVE_SET_GA_PROPERTY),Y)
      LIB_DEFINES += -DHAVE_SET_GA_PROPERTY
endif
