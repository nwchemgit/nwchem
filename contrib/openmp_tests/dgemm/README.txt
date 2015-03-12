README DGEMM
============

The aim of the current parallel DGEMM implementations is to compare issues
that may arise in different implementation approaches. In particular we are
interested in hybrid programming models <distributed>/OpenMP where <distributed>
is a distributed memory programming model that might be e.g. Global Arrays or 
MPI. These hybrid programming models are considered prime candidates for
programming shared memory cluster computers. The question is which paradigm 
fits the needs of NWChem as a highly scalable computational chemistry
application best?

In the kernels implemented we compute a matrix-matrix product C = A*B. The
matrices are square and partitioned onto the a square node grid. Within a
node the data may be further partitioned onto a square processor grid as well.
All of the blocks are also square and of equal dimensions. For convenience
we will typically use data locality of the blocks of matrix C to decide who
executes the associated work.

Initially instances A1 and B1 of the matrices A and B will be created on the
root processor. These matrices will be filled with random data and the
results distributed into matrices A and B. The parallel routine will compute
the distributed matrix C. Subsequently the root processor will compute C1
in serial and C1-C will be computed to check the correctness of the result.
The elapse time for the parallel dgemm will be measured and reported.

The list of tested implementations is:

dgemm_mpi.F        : An MPI only implementation
dgemm_ga.F         : A Global Array only implementation
dgemm_mpi_openmp.F : An MPI/OpenMP implementation
dgemm_ga_openmp.F  : A Global Array/OpenMP implementation

