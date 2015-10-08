# $Id$

         LIBRARY = libnwxc.a

    OBJ_OPTIMIZE = nwxc_print.o nwxc.o nwxc_add.o nwxc_block.o nwxc_query.o \
		   nwxc_nwchem.o nwxcm_eval.o \
                   nwxc_vdw_main.o nwxc_vdw_e.o nwxc_vdw_der.o nwxc_vdw_hess.o \
                   nwxc_vdw_util.o nwxc_vdw3a.o nwxc_vdw3b.o nwxc_vdw3c.o \
                   nwxcm_heaviside.o
                   #nwxc_fnc_print.o nwxcm_heaviside.o
#   OBJ_OPTIMIZE = nwad.o nwxc_eval.o nwxc_x_att.o nwxc_x_b88.o nwxc_x_bnl.o \
#                  nwxc_x_camb88.o nwxc_x_camlsd.o nwxc_x_campbe.o nwxc_x_camrpbe.o nwxc_x_cams12.o \
#                  nwxc_x_dirac.o nwxc_x_gill.o nwxc_x_hse08.o nwxc_x_ft97.o nwxc_x_m05.o nwxc_x_m06.o \
#                  nwxc_x_m11.o nwxc_x_opt.o nwxc_x_pbe.o nwxc_x_rpbe.o \
#                  nwxc_x_pkzb99.o nwxc_x_s12.o nwxc_x_sogga.o nwxc_x_tpss03.o \
#                  nwxc_x_wpbe.o nwxc_x_vs98.o nwxc_c_b95.o nwxc_c_ft97.o nwxc_c_lyp.o \
#                  nwxc_c_m05.o nwxc_c_m06.o nwxc_c_m11.o nwxc_c_opt.o \
#                  nwxc_c_p91.o nwxc_c_p91_vwn5.o nwxc_c_pw91lda.o nwxc_c_sogga.o \
#                  nwxc_c_pbe.o nwxc_c_Mpbe.o nwxc_c_pkzb99.o nwxc_c_tpss03.o \
#                  nwxc_c_vs98.o nwxc_c_lsda.o nwxc_c_vwn1.o \
#                  nwxc_c_vwn1_rpa.o nwxc_c_vwn2.o nwxc_c_vwn3.o nwxc_c_vwn4.o \
#                  nwxc_c_vwn5.o nwxc_c_op.o \
#                  nwxc_vdw_main.o nwxc_vdw_e.o nwxc_vdw_der.o nwxc_vdw_hess.o \
#                  nwxc_vdw_util.o nwxc_cs00.o nwxc_vdw3a.o nwxc_vdw3b.o nwxc_vdw3c.o \
#	   nwxcP_xc_b97.o nwxc_x_b97.o nwxc_c_b97.o nwxc_xc_kt1.o \
#                  nwxc_c_spbe.o nwxc_x_ssbD_1.o nwxc_c_perdew81.o \
#	   nwxc_c_perdew86.o nwxc_x_pw91.o nwxc_x_pw6.o nwxc_x_b86b.o 

#            OBJ = nwxc_block.o nwxc_add.o nwxc_query.o nwxc_nwchem.o \
#                  nwxc_vdw_init.o nwxc_vdw3d.o
             OBJ = nwxc_vdw_init.o nwxc_vdw3d.o

         HEADERS = nwxc.fh nwxc_param.fh

         SUBDIRS = maxima

    LIB_INCLUDES = -I./nwad

     LIB_DEFINES += -DNODEBUG
ifdef NWCHEM_NWXC_QUAD_PREC
     LIB_DEFINES += -DNWAD_QUAD_PREC -DNWXC_QUAD_PREC
endif
ifdef USE_FORTRAN2008
     LIB_DEFINES += -DUSE_FORTRAN2008
endif

     LIB_TARGETS = 

       TEST_LIBS = $(LIBRARY) $(LIBS) 

       USES_BLAS =  nwxc_vdw_main.F nwxc_vdw_e.F nwxc_vdw_der.F nwxc_vdw_hess.F nwxc_vdw_util.F

include ../config/makefile.h
include ../config/makelib.h

# NWXC_FNC_PRINT
# ==============
#
# The nwxc_fnc_print target is specifically meant to print density
# functionals as Maxima expressions. This target is not meant to do
# any real computation. Therefore it should not be built as part of
# routine compilations, but only upon a specific request.
#
# Also the are some special requirements to build this target. 
# 1. It builds only with the internal BLAS libraries.
# 2. You must build NWChem the regular first before attempting to 
#    build this target. (There are many dependencies on include files
#    and libraries).
#
# So the recommended approach to generating this target is to follow
# the recipe:
# 1. export BLASOPT=" "
# 2. run ./contrib/distro-tools/build_nwchem realclean
# 3. run ./contrib/distro-tools/build_nwchem
# 4. source ./nwchem_make_env.sh
# 5. cd src/nwxc
# 6. make nwxc_fnc_print
# 
.PRECIOUS: nwxc_fnc_print

nwxc_fnc_print: nwxc_fnc_print.F nwxc_eval_print.F stubs.F ../../lib/$(NWCHEM_TARGET)/$(LIBRARY)
	@(env | sort)
	$(FC) $(FFLAGS) $(INCLUDES) $(DEFINES) -o nwxc_fnc_print nwxc_fnc_print.F nwxc_eval_print.F stubs.F $(LDFLAGS) -lnwctask -lnwxc -lcons -lnwcutil -lga -larmci -lperfm -lcons -lnwcutil -llapack -lblas -L$(MPI_LIB) $(LIBMPI) $(COMM_LIBS)
#	$(LINK.f) -o nwxc_fnc_print nwxc_fnc_print.o -lnwctask -lnwxc -lcons -lnwcutil -lga -larmci -lperfm -lcons -lnwcutil -llapack -lblas $(MPI_LIB) $(LIBMPI) 

#nwad.fh: nwad0.mod nwad1.mod nwad2.mod nwad3.mod
#	touch nwad.fh

#nwad0.mod nwad1.mod nwad2.mod nwad3.mod nwad.o: ./nwad/nwad.F90
#	$(FC) -c $(FFLAGS) ./nwad/nwad.F90

#nwad.fh: ./nwad/nwad.F90
#	$(FC) -c $(FFLAGS) ./nwad/nwad.F90
#	touch nwad.fh

nwxc.o:	nwxc_print.o
