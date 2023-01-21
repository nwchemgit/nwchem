#
# $Id$
#
    OBJ_OPTIMIZE = nwpwxc_eval.o nwpwxc_x_att.o nwpwxc_x_b88.o nwpwxc_x_bnl.o \
                   nwpwxc_x_camb88.o nwpwxc_x_camlsd.o nwpwxc_x_campbe.o nwpwxc_x_cams12.o \
                   nwpwxc_x_dirac.o nwpwxc_x_gill.o nwpwxc_x_hse08.o nwpwxc_x_ft97.o nwpwxc_x_m05.o nwpwxc_x_m06.o \
                   nwpwxc_x_m11.o nwpwxc_x_opt.o nwpwxc_x_pbe.o nwpwxc_x_rpbe.o \
                   nwpwxc_x_pkzb99.o nwpwxc_x_s12.o nwpwxc_x_sogga.o nwpwxc_x_tpss03.o \
                   nwpwxc_x_wpbe.o nwpwxc_x_vs98.o nwpwxc_c_b95.o nwpwxc_c_ft97.o nwpwxc_c_lyp.o \
                   nwpwxc_c_m05.o nwpwxc_c_m06.o nwpwxc_c_m11.o nwpwxc_c_opt.o \
                   nwpwxc_c_p91.o nwpwxc_c_p91_vwn5.o nwpwxc_c_pw91lda.o nwpwxc_c_sogga.o \
                   nwpwxc_c_pbe.o nwpwxc_c_Mpbe.o nwpwxc_c_pkzb99.o nwpwxc_c_tpss03.o \
                   nwpwxc_c_vs98.o nwpwxc_c_lsda.o nwpwxc_c_vwn1.o \
                   nwpwxc_c_vwn1_rpa.o nwpwxc_c_vwn2.o nwpwxc_c_vwn3.o nwpwxc_c_vwn4.o \
                   nwpwxc_c_vwn5.o nwpwxc_c_op.o \
                   nwpwxc_cs00.o nwpwxc_vdw3b.o nwpwxc_vdw3c.o \
                   nwpwxcP_xc_b97.o nwpwxc_x_b97.o nwpwxc_c_b97.o nwpwxc_xc_kt1.o \
                   nwpwxc_c_spbe.o nwpwxc_x_ssbD_1.o nwpwxc_c_perdew81.o \
                   nwpwxc_c_perdew86.o nwpwxc_x_pw91.o nwpwxc_x_pw6.o nwpwxc_x_b86b.o \
                   nwpwxc_x_scan.o nwpwxc_c_scan.o
ifdef NO_NWPWXC_VDW3A
             OBJ_OPTIMIZE += nwpwxc_vdw3a_stub.o
else
             OBJ_OPTIMIZE += nwpwxc_vdw3a.o
endif

             OBJ = nwpwxc_block.o nwpwxc_add.o nwpwxc_query.o nwpwxc_nwchem.o \
                   nwpwxc_vdw3d.o

      LIBRARY = libnwpwlib.a

    USES_BLAS = nwpwxc_vdw3c.F nwpwxc_eval.F nwpwxc_c_vwn_a.fh
ifndef NO_NWPWXC_VDW3A
    USES_BLAS += nwpwxc_vdw3a.F 
endif

      HEADERS = nwpwxc.fh nwpwxc_param.fh

  LIB_DEFINES = 

include ../../../config/makefile.h
include ../../../config/makelib.h
ifndef NO_NWPWXC_VDW3A
nwpwxc_vdw3a.F:
	./build_dftd3a.sh
endif


