all: test_sin.log test_cos.log test_tan.log test_asin.log test_acos.log test_atan.log

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

nwad.mod: ../nwad.F90
	gfortran -c $^

%.o:	%.F
	gfortran -c $^

%.F:	%.max
	maxima -b $^
