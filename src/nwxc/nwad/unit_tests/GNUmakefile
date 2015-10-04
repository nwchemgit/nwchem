# $Id$
#
# Unit Tests for the NWAD Automatic Differentiation code
# ======================================================
#
# Running these unit tests requires the Maxima symbolic algebra package.
# See http://maxima.sourceforge.net/ for details. Maxima is used to 
# derive the reference implementations of the tests.
#
# Every test case, when run, produces a log-file listing all the errors.
# Hence empty log-files are a good thing. The labels on log-lines indicate
# which level of derivative shows problems:
# - F   : the basic expression is evaluated incorrectly
# - DF  : the 1st derivative expression is evaluated incorrectly
# - DF2 : the 2nd derivative expression is evaluated incorrectly
# - DF3 : the 3rd derivative expression is evaluated incorrectly
# The code compares against relative errors.
#
.PRECIOUS:
.PRECIOUS: .F

all: test_sin_m3.log \
     test_abs_m3.log \
     test_add_m3.log  test_addx_m3.log test_addy_m3.log \
     test_acos_m3.log test_asin_m3.log test_asinh_m3.log test_atan_m3.log\
     test_cos_m3.log  test_cosh_m3.log \
     test_div_m3.log  test_divx_m3.log test_divy_m3.log \
     test_erf_m3.log  test_erfc_m3.log test_exp_m3.log \
     test_mul_m3.log  test_mulx_m3.log test_muly_m3.log \
     test_log_m3.log \
     test_pow_m3.log  test_powx_m3.log test_powy_m3.log \
     test_sinh_m3.log test_sqrt_m3.log \
     test_sub_m3.log  test_subx_m3.log test_suby_m3.log \
     test_tan_m3.log  test_tanh_m3.log \
     test_sin_p3.log  test_cos_p3.log  test_tan_p3.log

#all: test_sin.log   test_cos.log  test_tan.log \
#     test_asin.log  test_acos.log test_atan.log \
#     test_sinh.log  test_cosh.log test_tanh.log \
#     test_asinh.log test_erf.log  test_erfc.log \
#     test_abs.log   test_exp.log  test_sqrt.log  \
#     test_log.log \
#     test_add.log   test_addx.log test_addy.log \
#     test_sub.log   test_subx.log test_suby.log \
#     test_mul.log   test_mulx.log test_muly.log \
#     test_div.log   test_divx.log test_divy.log \
#     test_pow.log   test_powx.log test_powy.log \
#     test_dx.log    test_dxy.log  test_dxyz.log \
#     test_lt.log    test_le.log   test_gt.log   test_ge.log

test_sin.log: test_sin
	echo "Running " $^
	./test_sin | tee test_sin.log

test_sin_m3.log: test_sin_m3
	echo "Running " $^
	./test_sin_m3 | tee test_sin_m3.log

test_sin_p3.log: test_sin_p3
	echo "Running " $^
	./test_sin_p3 | tee test_sin_p3.log

test_cos.log: test_cos
	echo "Running " $^
	./test_cos | tee test_cos.log

test_cos_m3.log: test_cos_m3
	echo "Running " $^
	./test_cos_m3 | tee test_cos_m3.log

test_cos_p3.log: test_cos_p3
	echo "Running " $^
	./test_cos_p3 | tee test_cos_p3.log

test_tan.log: test_tan
	echo "Running " $^
	./test_tan | tee test_tan.log

test_tan_m3.log: test_tan_m3
	echo "Running " $^
	./test_tan_m3 | tee test_tan_m3.log

test_tan_p3.log: test_tan_p3
	echo "Running " $^
	./test_tan_p3 | tee test_tan_p3.log

test_asin.log: test_asin
	echo "Running " $^
	./test_asin | tee test_asin.log

test_asin_m3.log: test_asin_m3
	echo "Running " $^
	./test_asin_m3 | tee test_asin_m3.log

test_acos.log: test_acos
	echo "Running " $^
	./test_acos | tee test_acos.log

test_acos_m3.log: test_acos_m3
	echo "Running " $^
	./test_acos_m3 | tee test_acos_m3.log

test_atan.log: test_atan
	echo "Running " $^
	./test_atan | tee test_atan.log

test_atan_m3.log: test_atan_m3
	echo "Running " $^
	./test_atan_m3 | tee test_atan_m3.log

test_sinh.log: test_sinh
	echo "Running " $^
	./test_sinh | tee test_sinh.log

test_asinh_m3.log: test_asinh_m3
	echo "Running " $^
	./test_asinh_m3 | tee test_asinh_m3.log

test_sinh_m3.log: test_sinh_m3
	echo "Running " $^
	./test_sinh_m3 | tee test_sinh_m3.log

test_cosh.log: test_cosh
	echo "Running " $^
	./test_cosh | tee test_cosh.log

test_cosh_m3.log: test_cosh_m3
	echo "Running " $^
	./test_cosh_m3 | tee test_cosh_m3.log

test_tanh.log: test_tanh
	echo "Running " $^
	./test_tanh | tee test_tanh.log

test_tanh_m3.log: test_tanh_m3
	echo "Running " $^
	./test_tanh_m3 | tee test_tanh_m3.log

test_asinh.log: test_asinh
	echo "Running " $^
	./test_asinh | tee test_asinh.log

test_erf.log: test_erf
	echo "Running " $^
	./test_erf | tee test_erf.log

test_erf_m3.log: test_erf_m3
	echo "Running " $^
	./test_erf_m3 | tee test_erf_m3.log

test_erfc.log: test_erfc
	echo "Running " $^
	./test_erfc | tee test_erfc.log

test_erfc_m3.log: test_erfc_m3
	echo "Running " $^
	./test_erfc_m3 | tee test_erfc_m3.log

test_abs.log: test_abs
	echo "Running " $^
	./test_abs | tee test_abs.log

test_abs_m3.log: test_abs_m3
	echo "Running " $^
	./test_abs_m3 | tee test_abs_m3.log

test_exp.log: test_exp
	echo "Running " $^
	./test_exp | tee test_exp.log

test_exp_m3.log: test_exp_m3
	echo "Running " $^
	./test_exp_m3 | tee test_exp_m3.log

test_sqrt.log: test_sqrt
	echo "Running " $^
	./test_sqrt | tee test_sqrt.log

test_sqrt_m3.log: test_sqrt_m3
	echo "Running " $^
	./test_sqrt_m3 | tee test_sqrt_m3.log

test_log.log: test_log
	echo "Running " $^
	./test_log | tee test_log.log

test_log_m3.log: test_log_m3
	echo "Running " $^
	./test_log_m3 | tee test_log_m3.log

test_add.log: test_add
	echo "Running " $^
	./test_add | tee test_add.log

test_add_m3.log: test_add_m3
	echo "Running " $^
	./test_add_m3 | tee test_add_m3.log

test_addx.log: test_addx
	echo "Running " $^
	./test_addx | tee test_addx.log

test_addx_m3.log: test_addx_m3
	echo "Running " $^
	./test_addx_m3 | tee test_addx_m3.log

test_addy.log: test_addy
	echo "Running " $^
	./test_addy | tee test_addy.log

test_addy_m3.log: test_addy_m3
	echo "Running " $^
	./test_addy_m3 | tee test_addy_m3.log

test_sub.log: test_sub
	echo "Running " $^
	./test_sub | tee test_sub.log

test_sub_m3.log: test_sub_m3
	echo "Running " $^
	./test_sub_m3 | tee test_sub_m3.log

test_subx.log: test_subx
	echo "Running " $^
	./test_subx | tee test_subx.log

test_subx_m3.log: test_subx_m3
	echo "Running " $^
	./test_subx_m3 | tee test_subx_m3.log

test_suby.log: test_suby
	echo "Running " $^
	./test_suby | tee test_suby.log

test_suby_m3.log: test_suby_m3
	echo "Running " $^
	./test_suby_m3 | tee test_suby_m3.log

test_mul.log: test_mul
	echo "Running " $^
	./test_mul | tee test_mul.log

test_mul_m3.log: test_mul_m3
	echo "Running " $^
	./test_mul_m3 | tee test_mul_m3.log

test_mulx.log: test_mulx
	echo "Running " $^
	./test_mulx | tee test_mulx.log

test_mulx_m3.log: test_mulx_m3
	echo "Running " $^
	./test_mulx_m3 | tee test_mulx_m3.log

test_muly.log: test_muly
	echo "Running " $^
	./test_muly | tee test_muly.log

test_muly_m3.log: test_muly_m3
	echo "Running " $^
	./test_muly_m3 | tee test_muly_m3.log

test_div.log: test_div
	echo "Running " $^
	./test_div | tee test_div.log

test_div_m3.log: test_div_m3
	echo "Running " $^
	./test_div_m3 | tee test_div_m3.log

test_divx.log: test_divx
	echo "Running " $^
	./test_divx | tee test_divx.log

test_divx_m3.log: test_divx_m3
	echo "Running " $^
	./test_divx_m3 | tee test_divx_m3.log

test_divy.log: test_divy
	echo "Running " $^
	./test_divy | tee test_divy.log

test_divy_m3.log: test_divy_m3
	echo "Running " $^
	./test_divy_m3 | tee test_divy_m3.log

test_pow.log: test_pow
	echo "Running " $^
	./test_pow | tee test_pow.log

test_pow_m3.log: test_pow_m3
	echo "Running " $^
	./test_pow_m3 | tee test_pow_m3.log

test_powx.log: test_powx
	echo "Running " $^
	./test_powx | tee test_powx.log

test_powx_m3.log: test_powx_m3
	echo "Running " $^
	./test_powx_m3 | tee test_powx_m3.log

test_powy.log: test_powy
	echo "Running " $^
	./test_powy | tee test_powy.log

test_powy_m3.log: test_powy_m3
	echo "Running " $^
	./test_powy_m3 | tee test_powy_m3.log

test_dx.log: test_dx
	echo "Running " $^
	./test_dx | tee test_dx.log

test_dxy.log: test_dxy
	echo "Running " $^
	./test_dxy | tee test_dxy.log

test_dxyz.log: test_dxyz
	echo "Running " $^
	./test_dxyz | tee test_dxyz.log

test_lt.log: test_lt
	echo "Running " $^
	./test_lt | tee test_lt.log

test_le.log: test_le
	echo "Running " $^
	./test_le | tee test_le.log

test_gt.log: test_gt
	echo "Running " $^
	./test_gt | tee test_gt.log

test_ge.log: test_ge
	echo "Running " $^
	./test_ge | tee test_ge.log

test_sin: nwad3.mod submaxima_sin.o test_sin.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o
	$(FC) test_sin.o submaxima_sin.o nwad.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o -o test_sin

test_sin_m3: nwad3.mod submaxima_sin3.o test_sin_m3.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o
	$(FC) test_sin_m3.o submaxima_sin3.o nwad.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o -o test_sin_m3

test_sin_p3: nwad3.mod submaxima_sin_p3.o test_sin_p3.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o
	$(FC) test_sin_p3.o submaxima_sin_p3.o nwad.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o -o test_sin_p3

test_cos: nwad3.mod submaxima_cos.o test_cos.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o
	$(FC) test_cos.o submaxima_cos.o nwad.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o -o test_cos

test_cos_m3: nwad3.mod submaxima_cos3.o test_cos_m3.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o
	$(FC) test_cos_m3.o submaxima_cos3.o nwad.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o -o test_cos_m3

test_cos_p3: nwad3.mod submaxima_cos_p3.o test_cos_p3.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o
	$(FC) test_cos_p3.o submaxima_cos_p3.o nwad.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o -o test_cos_p3

test_tan: nwad3.mod submaxima_tan.o test_tan.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o
	$(FC) test_tan.o submaxima_tan.o nwad.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o -o test_tan

test_tan_m3: nwad3.mod submaxima_tan3.o test_tan_m3.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o
	$(FC) test_tan_m3.o submaxima_tan3.o nwad.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o -o test_tan_m3

test_tan_p3: nwad3.mod submaxima_tan_p3.o test_tan_p3.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o
	$(FC) test_tan_p3.o submaxima_tan_p3.o nwad.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o -o test_tan_p3

test_asin: nwad3.mod submaxima_asin.o test_asin.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o
	$(FC) test_asin.o submaxima_asin.o nwad.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o -o test_asin

test_asin_m3: nwad3.mod submaxima_asin3.o test_asin_m3.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o
	$(FC) test_asin_m3.o submaxima_asin3.o nwad.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o -o test_asin_m3

test_acos: nwad3.mod submaxima_acos.o test_acos.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o
	$(FC) test_acos.o submaxima_acos.o nwad.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o -o test_acos

test_acos_m3: nwad3.mod submaxima_acos3.o test_acos_m3.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o
	$(FC) test_acos_m3.o submaxima_acos3.o nwad.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o -o test_acos_m3

test_atan: nwad3.mod submaxima_atan.o test_atan.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o
	$(FC) test_atan.o submaxima_atan.o nwad.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o -o test_atan

test_atan_m3: nwad3.mod submaxima_atan3.o test_atan_m3.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o
	$(FC) test_atan_m3.o submaxima_atan3.o nwad.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o -o test_atan_m3

test_sinh: nwad3.mod submaxima_sinh.o test_sinh.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o
	$(FC) test_sinh.o submaxima_sinh.o nwad.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o -o test_sinh

test_sinh_m3: nwad3.mod submaxima_sinh3.o test_sinh_m3.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o
	$(FC) test_sinh_m3.o submaxima_sinh3.o nwad.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o -o test_sinh_m3

test_cosh: nwad3.mod submaxima_cosh.o test_cosh.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o
	$(FC) test_cosh.o submaxima_cosh.o nwad.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o -o test_cosh

test_cosh_m3: nwad3.mod submaxima_cosh3.o test_cosh_m3.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o
	$(FC) test_cosh_m3.o submaxima_cosh3.o nwad.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o -o test_cosh_m3

test_tanh: nwad3.mod submaxima_tanh.o test_tanh.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o
	$(FC) test_tanh.o submaxima_tanh.o nwad.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o -o test_tanh

test_tanh_m3: nwad3.mod submaxima_tanh3.o test_tanh_m3.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o
	$(FC) test_tanh_m3.o submaxima_tanh3.o nwad.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o -o test_tanh_m3

test_asinh: nwad3.mod submaxima_asinh.o test_asinh.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o
	$(FC) test_asinh.o submaxima_asinh.o nwad.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o -o test_asinh

test_asinh_m3: nwad3.mod submaxima_asinh3.o test_asinh_m3.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o
	$(FC) test_asinh_m3.o submaxima_asinh3.o nwad.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o -o test_asinh_m3

test_erf: nwad3.mod submaxima_erf.o test_erf.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o
	$(FC) test_erf.o submaxima_erf.o nwad.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o -o test_erf

test_erf_m3: nwad3.mod submaxima_erf3.o test_erf_m3.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o
	$(FC) test_erf_m3.o submaxima_erf3.o nwad.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o -o test_erf_m3

test_erfc: nwad3.mod submaxima_erfc.o test_erfc.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o
	$(FC) test_erfc.o submaxima_erfc.o nwad.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o -o test_erfc

test_erfc_m3: nwad3.mod submaxima_erfc3.o test_erfc_m3.o errquit.o util_erf.o util_intrsc2008.o util_erf.o inp_strlen.o
	$(FC) test_erfc_m3.o submaxima_erfc3.o nwad.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o -o test_erfc_m3

test_abs: nwad3.mod submaxima_abs.o test_abs.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o
	$(FC) test_abs.o submaxima_abs.o nwad.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o -o test_abs

test_abs_m3: nwad3.mod submaxima_abs3.o test_abs_m3.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o
	$(FC) test_abs_m3.o submaxima_abs3.o nwad.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o -o test_abs_m3

test_exp: nwad3.mod submaxima_exp.o test_exp.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o
	$(FC) test_exp.o submaxima_exp.o nwad.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o -o test_exp

test_exp_m3: nwad3.mod submaxima_exp3.o test_exp_m3.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o
	$(FC) test_exp_m3.o submaxima_exp3.o nwad.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o -o test_exp_m3

test_sqrt: nwad3.mod submaxima_sqrt.o test_sqrt.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o
	$(FC) test_sqrt.o submaxima_sqrt.o nwad.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o -o test_sqrt

test_sqrt_m3: nwad3.mod submaxima_sqrt3.o test_sqrt_m3.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o
	$(FC) test_sqrt_m3.o submaxima_sqrt3.o nwad.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o -o test_sqrt_m3

test_log: nwad3.mod submaxima_log.o test_log.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o
	$(FC) test_log.o submaxima_log.o nwad.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o -o test_log

test_log_m3: nwad3.mod submaxima_log3.o test_log_m3.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o
	$(FC) test_log_m3.o submaxima_log3.o nwad.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o -o test_log_m3

test_add: nwad3.mod submaxima_add.o test_add.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o
	$(FC) test_add.o submaxima_add.o nwad.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o -o test_add

test_add_m3: nwad3.mod submaxima_add3.o test_add_m3.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o
	$(FC) test_add_m3.o submaxima_add3.o nwad.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o -o test_add_m3

test_addx: nwad3.mod submaxima_addx.o test_addx.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o
	$(FC) test_addx.o submaxima_addx.o nwad.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o -o test_addx

test_addx_m3: nwad3.mod submaxima_addx3.o test_addx_m3.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o
	$(FC) test_addx_m3.o submaxima_addx3.o nwad.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o -o test_addx_m3

test_addy: nwad3.mod submaxima_addy.o test_addy.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o
	$(FC) test_addy.o submaxima_addy.o nwad.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o -o test_addy

test_addy_m3: nwad3.mod submaxima_addy3.o test_addy_m3.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o
	$(FC) test_addy_m3.o submaxima_addy3.o nwad.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o -o test_addy_m3

test_sub: nwad3.mod submaxima_sub.o test_sub.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o
	$(FC) test_sub.o submaxima_sub.o nwad.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o -o test_sub

test_sub_m3: nwad3.mod submaxima_sub3.o test_sub_m3.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o
	$(FC) test_sub_m3.o submaxima_sub3.o nwad.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o -o test_sub_m3

test_subx: nwad3.mod submaxima_subx.o test_subx.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o
	$(FC) test_subx.o submaxima_subx.o nwad.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o -o test_subx

test_subx_m3: nwad3.mod submaxima_subx3.o test_subx_m3.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o
	$(FC) test_subx_m3.o submaxima_subx3.o nwad.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o -o test_subx_m3

test_suby: nwad3.mod submaxima_suby.o test_suby.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o
	$(FC) test_suby.o submaxima_suby.o nwad.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o -o test_suby

test_suby_m3: nwad3.mod submaxima_suby3.o test_suby_m3.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o
	$(FC) test_suby_m3.o submaxima_suby3.o nwad.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o -o test_suby_m3

test_mul: nwad3.mod submaxima_mul.o test_mul.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o
	$(FC) test_mul.o submaxima_mul.o nwad.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o -o test_mul

test_mul_m3: nwad3.mod submaxima_mul3.o test_mul_m3.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o
	$(FC) test_mul_m3.o submaxima_mul3.o nwad.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o -o test_mul_m3

test_mulx: nwad3.mod submaxima_mulx.o test_mulx.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o
	$(FC) test_mulx.o submaxima_mulx.o nwad.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o -o test_mulx

test_mulx_m3: nwad3.mod submaxima_mulx3.o test_mulx_m3.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o
	$(FC) test_mulx_m3.o submaxima_mulx3.o nwad.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o -o test_mulx_m3

test_muly: nwad3.mod submaxima_muly.o test_muly.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o
	$(FC) test_muly.o submaxima_muly.o nwad.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o -o test_muly

test_muly_m3: nwad3.mod submaxima_muly3.o test_muly_m3.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o
	$(FC) test_muly_m3.o submaxima_muly3.o nwad.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o -o test_muly_m3

test_div: nwad3.mod submaxima_div.o test_div.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o
	$(FC) test_div.o submaxima_div.o nwad.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o -o test_div

test_div_m3: nwad3.mod submaxima_div3.o test_div_m3.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o
	$(FC) test_div_m3.o submaxima_div3.o nwad.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o -o test_div_m3

test_divx: nwad3.mod submaxima_divx.o test_divx.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o
	$(FC) test_divx.o submaxima_divx.o nwad.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o -o test_divx

test_divx_m3: nwad3.mod submaxima_divx3.o test_divx_m3.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o
	$(FC) test_divx_m3.o submaxima_divx3.o nwad.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o -o test_divx_m3

test_divy: nwad3.mod submaxima_divy.o test_divy.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o
	$(FC) test_divy.o submaxima_divy.o nwad.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o -o test_divy

test_divy_m3: nwad3.mod submaxima_divy3.o test_divy_m3.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o
	$(FC) test_divy_m3.o submaxima_divy3.o nwad.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o -o test_divy_m3

test_pow: nwad3.mod submaxima_pow.o test_pow.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o
	$(FC) test_pow.o submaxima_pow.o nwad.o errquit.o -util_intrsc2008.o util_erf.o inp_strlen.o o test_pow

test_pow_m3: nwad3.mod submaxima_pow3.o test_pow_m3.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o
	$(FC) test_pow_m3.o submaxima_pow3.o nwad.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o -o test_pow_m3

test_powx: nwad3.mod submaxima_powx.o test_powx.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o
	$(FC) test_powx.o submaxima_powx.o nwad.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o -o test_powx

test_powx_m3: nwad3.mod submaxima_powx3.o test_powx_m3.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o
	$(FC) test_powx_m3.o submaxima_powx3.o nwad.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o -o test_powx_m3

test_powy: nwad3.mod submaxima_powy.o test_powy.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o
	$(FC) test_powy.o submaxima_powy.o nwad.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o -o test_powy

test_powy_m3: nwad3.mod submaxima_powy3.o test_powy_m3.o errquit.o
	$(FC) test_powy_m3.o submaxima_powy3.o nwad.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o -o test_powy_m3

test_dx: nwad3.mod submaxima_dx.o test_dx.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o
	$(FC) test_dx.o submaxima_dx.o nwad.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o -o test_dx

test_dxy: nwad3.mod submaxima_dxy.o test_dxy.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o
	$(FC) test_dxy.o submaxima_dxy.o nwad.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o -o test_dxy

test_dxyz: nwad3.mod submaxima_dxyz.o test_dxyz.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o
	$(FC) test_dxyz.o submaxima_dxyz.o nwad.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o -o test_dxyz

test_lt: nwad3.mod test_lt.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o
	$(FC) test_lt.o nwad.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o -o test_lt

test_le: nwad3.mod test_le.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o
	$(FC) test_le.o nwad.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o -o test_le

test_gt: nwad3.mod test_gt.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o
	$(FC) test_gt.o nwad.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o -o test_gt

test_ge: nwad3.mod test_ge.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o
	$(FC) test_ge.o nwad.o errquit.o util_intrsc2008.o util_erf.o inp_strlen.o -o test_ge

util_intrsc2008.o: ../../../util/util_intrsc2008.F
	$(FC) -c ../../../util/util_intrsc2008.F

util_erf.o: ../../../util/util_erf.F
	$(FC) -c ../../../util/util_erf.F

opt_powx: opt_powx.F ../nwad.F
	$(FC) -DUSE_FORTRAN2008 -DNWAD_DEG=3 opt_powx.F -o opt_powx

# NWAD_DEG not used yet, but planned as a way to generate modules for different
# degrees of differentiation
nwad3.mod: ../nwad.F
	$(FC) -I../../../inp -c $^

%.o:	%.F
	$(FC) -c $^

%.F:	%.max
	maxima -b $^
	tar -zcf subs.tgz submaxima_*.F

clean:
	rm -f *.o *.log submaxima_*.F nwad?.mod
	find . -maxdepth 1 -perm -u=x -type f -exec rm {} \; -print
