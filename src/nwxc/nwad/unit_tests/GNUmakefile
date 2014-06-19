all: test_sin.log test_cos.log test_tan.log \
     test_asin.log test_acos.log test_atan.log

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

nwad.mod: ../nwad.F90
	gfortran -c $^

%.o:	%.F
	gfortran -c $^

%.F:	%.max
	maxima -b $^

clean:
	rm -f *.o *.log submaxima_*.F nwad.mod
	find . -maxdepth 1 -perm -u=x -type f -exec rm {} \; -print
