# $Id: GNUmakefile 21484 2011-11-09 21:22:15Z d3p852 $

         LIBRARY = libnwxc.a

    OBJ_OPTIMIZE = nwxc_vdw_main.o nwxc_vdw_e.o nwxc_vdw_der.o nwxc_vdw_hess.o nwxc_vdw_util.o

             OBJ = nwxc_block.o nwxc_add.o nwxc_query.o nwxc_nwchem.o \
                   nwxc_x_att.o nwxc_x_b88.o nwxc_x_bnl.o \
                   nwxc_x_camb88.o nwxc_x_camlsd.o nwxc_x_campbe.o \
                   nwxc_x_dirac.o nwxc_x_hse08.o \
                   nwxc_x_pbe.o nwxc_x_rpbe.o nwxc_x_wpbe.o nwxc_x_vs98.o \
                   nwxc_c_b95.o \
                   nwxc_c_pw91lda.o nwxc_c_pbe.o nwxc_c_vs98.o \
                   nwxc_c_lsda.o nwxc_eval.o \
                   nwxc_c_vwn1.o nwxc_c_vwn1_rpa.o nwxc_c_vwn2.o \
                   nwxc_c_vwn3.o nwxc_c_vwn4.o nwxc_c_vwn5.o \
                   nwxc_vdw_init.o

         HEADERS = nwxc.fh nwxc_param.fh

     LIB_TARGETS = 

       TEST_LIBS = $(LIBRARY) $(LIBS) 

USES_BLAS =  nwxc_vdw_main.F nwxc_vdw_e.F nwxc_vdw_der.F nwxc_vdw_hess.F nwxc_vdw_util.F

include ../config/makefile.h
include ../config/makelib.h
