#include <mpi.h>
#ifdef XLFLINUX
void utilc_mpicommrank_(
#else
void utilc_mpicommrank(
#endif
		       int *me)
{
  int myid;
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);
  *me= myid;
}
#ifdef XLFLINUX
void utilc_mpicommsize_(
#else
void utilc_mpicommsize(
#endif
		       int *size)
{
  int mysize;
  MPI_Comm_size(MPI_COMM_WORLD,&mysize);
  *size= mysize;
}
