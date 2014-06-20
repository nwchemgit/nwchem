all: test_sin.log test_cos.log test_tan.log \
     test_asin.log test_acos.log test_atan.log \
     test_sinh.log test_cosh.log test_tanh.log \
     test_asinh.log test_erf.log test_erfc.log \
     test_abs.log test_exp.log test_sqrt.log  \
     test_log.log \
     test_add.log test_addx.log test_addy.log \
     test_sub.log test_subx.log test_suby.log \
     test_mul.log test_mulx.log test_muly.log \
     test_div.log test_divx.log test_divy.log \
     test_pow.log test_powx.log test_powy.log

test_sin.log: test_sin
	echo "Running " $^
	./test_sin | tee test_sin.log

test_cos.log: test_cos
	echo "Running " $^
	./test_cos | tee test_cos.log

test_tan.log: test_tan
	echo "Running " $^
	./test_tan | tee test_tan.log

test_asin.log: test_asin
	echo "Running " $^
	./test_asin | tee test_asin.log

test_acos.log: test_acos
	echo "Running " $^
	./test_acos | tee test_acos.log

test_atan.log: test_atan
	echo "Running " $^
	./test_atan | tee test_atan.log

test_sinh.log: test_sinh
	echo "Running " $^
	./test_sinh | tee test_sinh.log

test_cosh.log: test_cosh
	echo "Running " $^
	./test_cosh | tee test_cosh.log

test_tanh.log: test_tanh
	echo "Running " $^
	./test_tanh | tee test_tanh.log

test_asinh.log: test_asinh
	echo "Running " $^
	./test_asinh | tee test_asinh.log

test_erf.log: test_erf
	echo "Running " $^
	./test_erf | tee test_erf.log

test_erfc.log: test_erfc
	echo "Running " $^
	./test_erfc | tee test_erfc.log

test_abs.log: test_abs
	echo "Running " $^
	./test_abs | tee test_abs.log

test_exp.log: test_exp
	echo "Running " $^
	./test_exp | tee test_exp.log

test_sqrt.log: test_sqrt
	echo "Running " $^
	./test_sqrt | tee test_sqrt.log

test_log.log: test_log
	echo "Running " $^
	./test_log | tee test_log.log

test_add.log: test_add
	echo "Running " $^
	./test_add | tee test_add.log

test_addx.log: test_addx
	echo "Running " $^
	./test_addx | tee test_addx.log

test_addy.log: test_addy
	echo "Running " $^
	./test_addy | tee test_addy.log

test_sub.log: test_sub
	echo "Running " $^
	./test_sub | tee test_sub.log

test_subx.log: test_subx
	echo "Running " $^
	./test_subx | tee test_subx.log

test_suby.log: test_suby
	echo "Running " $^
	./test_suby | tee test_suby.log

test_mul.log: test_mul
	echo "Running " $^
	./test_mul | tee test_mul.log

test_mulx.log: test_mulx
	echo "Running " $^
	./test_mulx | tee test_mulx.log

test_muly.log: test_muly
	echo "Running " $^
	./test_muly | tee test_muly.log

test_div.log: test_div
	echo "Running " $^
	./test_div | tee test_div.log

test_divx.log: test_divx
	echo "Running " $^
	./test_divx | tee test_divx.log

test_divy.log: test_divy
	echo "Running " $^
	./test_divy | tee test_divy.log

test_pow.log: test_pow
	echo "Running " $^
	./test_pow | tee test_pow.log

test_powx.log: test_powx
	echo "Running " $^
	./test_powx | tee test_powx.log

test_powy.log: test_powy
	echo "Running " $^
	./test_powy | tee test_powy.log

test_sin: nwad.mod submaxima_sin.o test_sin.o
	gfortran test_sin.o submaxima_sin.o nwad.o -o test_sin

test_cos: nwad.mod submaxima_cos.o test_cos.o
	gfortran test_cos.o submaxima_cos.o nwad.o -o test_cos

test_tan: nwad.mod submaxima_tan.o test_tan.o
	gfortran test_tan.o submaxima_tan.o nwad.o -o test_tan

test_asin: nwad.mod submaxima_asin.o test_asin.o
	gfortran test_asin.o submaxima_asin.o nwad.o -o test_asin

test_acos: nwad.mod submaxima_acos.o test_acos.o
	gfortran test_acos.o submaxima_acos.o nwad.o -o test_acos

test_atan: nwad.mod submaxima_atan.o test_atan.o
	gfortran test_atan.o submaxima_atan.o nwad.o -o test_atan

test_sinh: nwad.mod submaxima_sinh.o test_sinh.o
	gfortran test_sinh.o submaxima_sinh.o nwad.o -o test_sinh

test_cosh: nwad.mod submaxima_cosh.o test_cosh.o
	gfortran test_cosh.o submaxima_cosh.o nwad.o -o test_cosh

test_tanh: nwad.mod submaxima_tanh.o test_tanh.o
	gfortran test_tanh.o submaxima_tanh.o nwad.o -o test_tanh

test_asinh: nwad.mod submaxima_asinh.o test_asinh.o
	gfortran test_asinh.o submaxima_asinh.o nwad.o -o test_asinh

test_erf: nwad.mod submaxima_erf.o test_erf.o
	gfortran test_erf.o submaxima_erf.o nwad.o -o test_erf

test_erfc: nwad.mod submaxima_erfc.o test_erfc.o
	gfortran test_erfc.o submaxima_erfc.o nwad.o -o test_erfc

test_abs: nwad.mod submaxima_abs.o test_abs.o
	gfortran test_abs.o submaxima_abs.o nwad.o -o test_abs

test_exp: nwad.mod submaxima_exp.o test_exp.o
	gfortran test_exp.o submaxima_exp.o nwad.o -o test_exp

test_sqrt: nwad.mod submaxima_sqrt.o test_sqrt.o
	gfortran test_sqrt.o submaxima_sqrt.o nwad.o -o test_sqrt

test_log: nwad.mod submaxima_log.o test_log.o
	gfortran test_log.o submaxima_log.o nwad.o -o test_log

test_add: nwad.mod submaxima_add.o test_add.o
	gfortran test_add.o submaxima_add.o nwad.o -o test_add

test_addx: nwad.mod submaxima_addx.o test_addx.o
	gfortran test_addx.o submaxima_addx.o nwad.o -o test_addx

test_addy: nwad.mod submaxima_addy.o test_addy.o
	gfortran test_addy.o submaxima_addy.o nwad.o -o test_addy

test_sub: nwad.mod submaxima_sub.o test_sub.o
	gfortran test_sub.o submaxima_sub.o nwad.o -o test_sub

test_subx: nwad.mod submaxima_subx.o test_subx.o
	gfortran test_subx.o submaxima_subx.o nwad.o -o test_subx

test_suby: nwad.mod submaxima_suby.o test_suby.o
	gfortran test_suby.o submaxima_suby.o nwad.o -o test_suby

test_mul: nwad.mod submaxima_mul.o test_mul.o
	gfortran test_mul.o submaxima_mul.o nwad.o -o test_mul

test_mulx: nwad.mod submaxima_mulx.o test_mulx.o
	gfortran test_mulx.o submaxima_mulx.o nwad.o -o test_mulx

test_muly: nwad.mod submaxima_muly.o test_muly.o
	gfortran test_muly.o submaxima_muly.o nwad.o -o test_muly

test_div: nwad.mod submaxima_div.o test_div.o
	gfortran test_div.o submaxima_div.o nwad.o -o test_div

test_divx: nwad.mod submaxima_divx.o test_divx.o
	gfortran test_divx.o submaxima_divx.o nwad.o -o test_divx

test_divy: nwad.mod submaxima_divy.o test_divy.o
	gfortran test_divy.o submaxima_divy.o nwad.o -o test_divy

test_pow: nwad.mod submaxima_pow.o test_pow.o
	gfortran test_pow.o submaxima_pow.o nwad.o -o test_pow

test_powx: nwad.mod submaxima_powx.o test_powx.o
	gfortran test_powx.o submaxima_powx.o nwad.o -o test_powx

test_powy: nwad.mod submaxima_powy.o test_powy.o
	gfortran test_powy.o submaxima_powy.o nwad.o -o test_powy

nwad.mod: ../nwad.F90
	gfortran -c $^

%.o:	%.F
	gfortran -c $^

%.F:	%.max
	maxima -b $^

clean:
	rm -f *.o *.log submaxima_*.F nwad.mod
	find . -maxdepth 1 -perm -u=x -type f -exec rm {} \; -print
