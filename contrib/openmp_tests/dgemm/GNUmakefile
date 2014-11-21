FC = mpif90

dgemm_mpi: dgemm_mpi.F dgemm.F xerbla.F lsame.F
	$(FC) -g -O0 -o dgemm_mpi dgemm_mpi.F dgemm.F xerbla.F lsame.F
