FNCPRINT = ../nwxc_fnc_print
AUTOXC   = ./bin/autoxc

LIBRARY = libnwxc.a

OBJ =  \
	f77/nwxcm_c_b95.o \
	f77/nwxcm_c_b97.o \
	f77/nwxcm_c_dldf.o \
	f77/nwxcm_c_ft97.o \
	f77/nwxcm_c_hcth.o \
	f77/nwxcm_c_lyp.o \
	f77/nwxcm_c_m05.o \
	f77/nwxcm_c_m06.o \
	f77/nwxcm_c_m11.o \
	f77/nwxcm_c_op.o \
	f77/nwxcm_c_opt.o \
	f77/nwxcm_c_p91.o \
	f77/nwxcm_c_p91_vwn5.o \
	f77/nwxcm_c_pbe.o \
	f77/nwxcm_c_perdew81.o \
	f77/nwxcm_c_perdew86.o \
	f77/nwxcm_c_pkzb99.o \
	f77/nwxcm_c_pw91lda.o \
	f77/nwxcm_c_sogga.o \
	f77/nwxcm_c_spbe.o \
	f77/nwxcm_c_tpss03.o \
	f77/nwxcm_c_vs98.o \
	f77/nwxcm_c_vwn1.o \
	f77/nwxcm_c_vwn1_rpa.o \
	f77/nwxcm_c_vwn2.o \
	f77/nwxcm_c_vwn3.o \
	f77/nwxcm_c_vwn4.o \
	f77/nwxcm_c_vwn5.o \
	f77/nwxcm_x_b86b.o \
	f77/nwxcm_x_b88.o \
	f77/nwxcm_x_b97.o \
	f77/nwxcm_x_bnl.o \
	f77/nwxcm_x_camb88.o \
	f77/nwxcm_x_camlsd.o \
	f77/nwxcm_x_campbe.o \
	f77/nwxcm_x_camrpbe.o \
	f77/nwxcm_x_cams12.o \
	f77/nwxcm_xc_kt1.o \
	f77/nwxcm_x_dirac.o \
	f77/nwxcm_x_dldf.o \
	f77/nwxcm_x_ft97.o \
	f77/nwxcm_x_gill.o \
	f77/nwxcm_x_hcth.o \
	f77/nwxcm_x_m05.o \
	f77/nwxcm_x_m06_2x.o \
	f77/nwxcm_x_m06.o \
	f77/nwxcm_x_m08.o \
	f77/nwxcm_x_m11.o \
	f77/nwxcm_x_opt.o \
	f77/nwxcm_x_pbe.o \
	f77/nwxcm_x_pkzb99.o \
	f77/nwxcm_x_pw6.o \
	f77/nwxcm_x_pw91.o \
	f77/nwxcm_x_rpbe.o \
	f77/nwxcm_x_s12.o \
	f77/nwxcm_x_sogga.o \
	f77/nwxcm_x_ssbd.o \
	f77/nwxcm_x_tpss03.o \
	f77/nwxcm_x_vs98.o \
	f77/nwxcm_x_wpbe.o \

include ../../config/makefile.h
include ../../config/makelib.h

f77/%.o: f77/%.F
	$(FC) $(FFLAGS) -c $< -o $@

f77/%.F: max/%.max
	$(AUTOXC) -o $@ $<

max/%.max: input/%.nw
	$(FNCPRINT) $<
