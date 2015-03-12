FC = mpif90
GA_LIB = -L../../../src/tools/install/lib
GA_INC = -I../../../src/tools/install/include

all: dgemm_mpi dgemm_ga

dgemm_mpi: dgemm_mpi.F dgemm.F xerbla.F lsame.F
	$(FC) -g -O0 -o dgemm_mpi dgemm_mpi.F dgemm.F xerbla.F lsame.F

dgemm_ga: dgemm_ga.F dgemm.F xerbla.F lsame.F ga_stubs.F
	$(FC) -g -O0 -i8 $(GA_INC) -o dgemm_ga dgemm_ga.F dgemm.F xerbla.F lsame.F ga_stubs.F $(GA_LIB) -lga -larmci
