FNCPRINT  = ../nwxc_fnc_print
ifndef DO_MAXIMA
AUTOXC    = echo
AUTOXC_DS = echo
else
# 
# To run MAXIMA uncomment these lines
# -----------------------------------
#
# Note that generating all the functional source code took 62 hours on 2 cores
# hence do this only when you have to! For example when you want to add a
# new functional run 
#
#	touch f77/*.F
#
# to update the timestamps on the pre-existing routines to avoid regenerating
# those. That way only the subroutine for your new functional will be generated.
#
# In addition note that if you need to generate the Maxima expression itself
# as well you will need the nwxc_fnc_print executable. To generate that run
#
#	cd ..
#	make nwxc_fnc_print
#
# There is no dependency for generating this executable as that would trigger
# regenerating all the *.max and *.F files. In addition this process is only
# occassionally needed so the inconvience of having to do this manually should
# be insignificant (famous last words...).
#
AUTOXC    = ./bin/autoxc
AUTOXC_DS = ./bin/autoxc-Ds
endif

LIBRARY = libnwxc.a

OBJ =  \
	nwxcm_c_b95.o \
	nwxcm_c_b97.o \
	nwxcm_c_bop.o \
	nwxcm_c_dldf.o \
	nwxcm_c_ft97.o \
	nwxcm_c_hcth.o \
	nwxcm_c_lyp.o \
	nwxcm_c_m05.o \
	nwxcm_c_m06.o \
	nwxcm_c_m11.o \
	nwxcm_c_mpbe.o \
	nwxcm_c_opt.o \
	nwxcm_c_p91.o \
	nwxcm_c_p91_vwn5.o \
	nwxcm_c_pbe.o \
	nwxcm_c_pbeop.o \
	nwxcm_c_perdew81.o \
	nwxcm_c_perdew86.o \
	nwxcm_c_pkzb99.o \
	nwxcm_c_pw91lda.o \
	nwxcm_c_sogga.o \
	nwxcm_c_sop.o \
	nwxcm_c_spbe.o \
	nwxcm_c_tpss03.o \
	nwxcm_c_vs98.o \
	nwxcm_c_vwn1.o \
	nwxcm_c_vwn1_rpa.o \
	nwxcm_c_vwn2.o \
	nwxcm_c_vwn3.o \
	nwxcm_c_vwn4.o \
	nwxcm_c_vwn5.o \
	nwxcm_x_b86b.o \
	nwxcm_x_b88.o \
	nwxcm_x_b97.o \
	nwxcm_x_bnl.o \
	nwxcm_x_camb88.o \
	nwxcm_x_camlsd.o \
	nwxcm_x_campbe.o \
	nwxcm_x_camrpbe.o \
	nwxcm_x_cams12.o \
	nwxcm_xc_kt1.o \
	nwxcm_x_dirac.o \
	nwxcm_x_dldf.o \
	nwxcm_x_ft97.o \
	nwxcm_x_gill.o \
	nwxcm_x_hcth.o \
	nwxcm_x_m05.o \
	nwxcm_x_m06_2x.o \
	nwxcm_x_m06.o \
	nwxcm_x_m08.o \
	nwxcm_x_m11.o \
	nwxcm_x_opt.o \
	nwxcm_x_pbe.o \
	nwxcm_x_pkzb99.o \
	nwxcm_x_pw6.o \
	nwxcm_x_pw91.o \
	nwxcm_x_rpbe.o \
	nwxcm_x_s12.o \
	nwxcm_x_sogga.o \
	nwxcm_x_ssbd.o \
	nwxcm_x_tpss03.o \
	nwxcm_x_vs98.o \
	nwxcm_x_wpbe.o


ifdef NWCHEM_NWXC_QUAD_PREC
     LIB_DEFINES = -DNWAD_QUAD_PREC -DNWXC_QUAD_PREC
endif

include ../../config/makefile.h
include ../../config/makelib.h

ifeq ($(_FC),gfortran)
  FDEBUG +=  -fno-gcse  -fno-tree-dominator-opts
  FDEBUG +=  --param ggc-min-expand=20 --param ggc-min-heapsize=65536
endif


# This Makefile should not be executed in parallel as the compilers
# will run out of memory in that case.
#
ifndef  USE_PARALLELMAKE
.NOTPARALLEL:
endif

.SECONDARY:
.SECONDARY: \
	max/nwxcm_c_b95.max \
	max/nwxcm_c_b97.max \
	max/nwxcm_c_bop.max \
	max/nwxcm_c_dldf.max \
	max/nwxcm_c_ft97.max \
	max/nwxcm_c_hcth.max \
	max/nwxcm_c_lyp.max \
	max/nwxcm_c_m05.max \
	max/nwxcm_c_m06.max \
	max/nwxcm_c_m11.max \
	max/nwxcm_c_mpbe.max \
	max/nwxcm_c_opt.max \
	max/nwxcm_c_p91.max \
	max/nwxcm_c_p91_vwn5.max \
	max/nwxcm_c_pbe.max \
	max/nwxcm_c_pbeop.max \
	max/nwxcm_c_perdew81.max \
	max/nwxcm_c_perdew86.max \
	max/nwxcm_c_pkzb99.max \
	max/nwxcm_c_pw91lda.max \
	max/nwxcm_c_sogga.max \
	max/nwxcm_c_sop.max \
	max/nwxcm_c_spbe.max \
	max/nwxcm_c_tpss03.max \
	max/nwxcm_c_vs98.max \
	max/nwxcm_c_vwn1.max \
	max/nwxcm_c_vwn1_rpa.max \
	max/nwxcm_c_vwn2.max \
	max/nwxcm_c_vwn3.max \
	max/nwxcm_c_vwn4.max \
	max/nwxcm_c_vwn5.max \
	max/nwxcm_x_b86b.max \
	max/nwxcm_x_b88.max \
	max/nwxcm_x_b97.max \
	max/nwxcm_x_bnl.max \
	max/nwxcm_x_camb88.max \
	max/nwxcm_x_camlsd.max \
	max/nwxcm_x_campbe.max \
	max/nwxcm_x_camrpbe.max \
	max/nwxcm_x_cams12.max \
	max/nwxcm_xc_kt1.max \
	max/nwxcm_x_dirac.max \
	max/nwxcm_x_dldf.max \
	max/nwxcm_x_ft97.max \
	max/nwxcm_x_gill.max \
	max/nwxcm_x_hcth.max \
	max/nwxcm_x_m05.max \
	max/nwxcm_x_m06_2x.max \
	max/nwxcm_x_m06.max \
	max/nwxcm_x_m08.max \
	max/nwxcm_x_m11.max \
	max/nwxcm_x_opt.max \
	max/nwxcm_x_pbe.max \
	max/nwxcm_x_pkzb99.max \
	max/nwxcm_x_pw6.max \
	max/nwxcm_x_pw91.max \
	max/nwxcm_x_rpbe.max \
	max/nwxcm_x_s12.max \
	max/nwxcm_x_sogga.max \
	max/nwxcm_x_ssbd.max \
	max/nwxcm_x_tpss03.max \
	max/nwxcm_x_vs98.max \
	max/nwxcm_x_wpbe.max \
	f77/nwxcm_c_b95.F \
	f77/nwxcm_c_b97.F \
	f77/nwxcm_c_bop.F \
	f77/nwxcm_c_dldf.F \
	f77/nwxcm_c_ft97.F \
	f77/nwxcm_c_hcth.F \
	f77/nwxcm_c_lyp.F \
	f77/nwxcm_c_m05.F \
	f77/nwxcm_c_m06.F \
	f77/nwxcm_c_m11.F \
	f77/nwxcm_c_mpbe.F \
	f77/nwxcm_c_opt.F \
	f77/nwxcm_c_p91.F \
	f77/nwxcm_c_p91_vwn5.F \
	f77/nwxcm_c_pbe.F \
	f77/nwxcm_c_pbeop.F \
	f77/nwxcm_c_perdew81.F \
	f77/nwxcm_c_perdew86.F \
	f77/nwxcm_c_pkzb99.F \
	f77/nwxcm_c_pw91lda.F \
	f77/nwxcm_c_sogga.F \
	f77/nwxcm_c_sop.F \
	f77/nwxcm_c_spbe.F \
	f77/nwxcm_c_tpss03.F \
	f77/nwxcm_c_vs98.F \
	f77/nwxcm_c_vwn1.F \
	f77/nwxcm_c_vwn1_rpa.F \
	f77/nwxcm_c_vwn2.F \
	f77/nwxcm_c_vwn3.F \
	f77/nwxcm_c_vwn4.F \
	f77/nwxcm_c_vwn5.F \
	f77/nwxcm_x_b86b.F \
	f77/nwxcm_x_b88.F \
	f77/nwxcm_x_b97.F \
	f77/nwxcm_x_bnl.F \
	f77/nwxcm_x_camb88.F \
	f77/nwxcm_x_camlsd.F \
	f77/nwxcm_x_campbe.F \
	f77/nwxcm_x_camrpbe.F \
	f77/nwxcm_x_cams12.F \
	f77/nwxcm_xc_kt1.F \
	f77/nwxcm_x_dirac.F \
	f77/nwxcm_x_dldf.F \
	f77/nwxcm_x_ft97.F \
	f77/nwxcm_x_gill.F \
	f77/nwxcm_x_hcth.F \
	f77/nwxcm_x_m05.F \
	f77/nwxcm_x_m06_2x.F \
	f77/nwxcm_x_m06.F \
	f77/nwxcm_x_m08.F \
	f77/nwxcm_x_m11.F \
	f77/nwxcm_x_opt.F \
	f77/nwxcm_x_pbe.F \
	f77/nwxcm_x_pkzb99.F \
	f77/nwxcm_x_pw6.F \
	f77/nwxcm_x_pw91.F \
	f77/nwxcm_x_rpbe.F \
	f77/nwxcm_x_s12.F \
	f77/nwxcm_x_sogga.F \
	f77/nwxcm_x_ssbd.F \
	f77/nwxcm_x_tpss03.F \
	f77/nwxcm_x_vs98.F \
	f77/nwxcm_x_wpbe.F \
	f77/nwxcm_c_vs98.F.start \
	f77/nwxcm_c_vs98.F.finished \
	f77/nwxcm_c_dldf.F.start \
	f77/nwxcm_c_dldf.F.finished \
	f77/nwxcm_c_m05.F.start \
	f77/nwxcm_c_m05.F.finished \
	f77/nwxcm_c_m06.F.start \
	f77/nwxcm_c_m06.F.finished \
	f77/nwxcm_c_mpbe.F.start \
	f77/nwxcm_c_mpbe.F.finished \
	f77/nwxcm_c_pbe.F.start \
	f77/nwxcm_c_pbe.F.finished \
	f77/nwxcm_c_pkzb99.F.start \
	f77/nwxcm_c_pkzb99.F.finished \
	f77/nwxcm_c_pw91lda.F.start \
	f77/nwxcm_c_pw91lda.F.finished \
	f77/nwxcm_c_tpss03.F.start \
	f77/nwxcm_c_tpss03.F.finished \
	f77/nwxcm_x_tpss03.F.start \
	f77/nwxcm_x_tpss03.F.finished



f77/nwxcm_c_vs98.F: max/nwxcm_c_vs98.max
	date > $@.start
	$(AUTOXC_DS) -o $@ $<
	date > $@.finished

f77/nwxcm_c_dldf.F: max/nwxcm_c_dldf.max
	date > $@.start
	$(AUTOXC_DS) -o $@ $<
	date > $@.finished

f77/nwxcm_c_m05.F: max/nwxcm_c_m05.max
	date > $@.start
	$(AUTOXC_DS) -o $@ $<
	date > $@.finished

f77/nwxcm_c_m06.F: max/nwxcm_c_m06.max
	date > $@.start
	$(AUTOXC_DS) -o $@ $<
	date > $@.finished

f77/nwxcm_c_mpbe.F: max/nwxcm_c_mpbe.max
	date > $@.start
	$(AUTOXC) -o $@ $<
	date > $@.finished

f77/nwxcm_c_pbe.F: max/nwxcm_c_pbe.max
	date > $@.start
	$(AUTOXC) -o $@ $<
	date > $@.finished

f77/nwxcm_c_pkzb99.F: max/nwxcm_c_pkzb99.max
	date > $@.start
	$(AUTOXC_DS) -o $@ $<
	date > $@.finished

f77/nwxcm_c_tpss03.F: max/nwxcm_c_tpss03.max
	date > $@.start
	$(AUTOXC_DS) -o $@ $<
	date > $@.finished

f77/nwxcm_x_tpss03.F: max/nwxcm_x_tpss03.max
	date > $@.start
	$(AUTOXC_DS) -o $@ $<
	date > $@.finished

%.o: f77/%.F
	$(FC) $(FFLAGS) $(LIB_DEFINES) -I./f77 -c $< 

f77/%.F: max/%.max
	$(AUTOXC) -o $@ $<

max/%.max: input/%.nw
	$(FNCPRINT) $<
