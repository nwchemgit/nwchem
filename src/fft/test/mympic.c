/*
 $Id: mympic.c,v 1.2 1997-11-04 10:07:55 d3e129 Exp $
 */

#include "mpi.h"


mysyncc()
{
int istatus;
MPI_Comm handle;
handle = MPI_COMM_WORLD;
istatus = mpi_barrier(handle);
}
