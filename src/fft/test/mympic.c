/*
 $Id$
 */

#include "mpi.h"


mysyncc()
{
int istatus;
MPI_Comm handle;
handle = MPI_COMM_WORLD;
istatus = mpi_barrier(handle);
}
