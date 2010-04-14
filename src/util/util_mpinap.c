/*
 $Id$
 */
#include <stdio.h>
#ifdef MPI
#include <mpi.h>
#else
#  include <global.h>
#endif
#include "global.h"
#include "macdecls.h"
/* This makes each process sleep for Mpirank/factor microsec*/
/*#define DEBUG 1*/
void FATR util_mpinap_(Integer *factor)
{
  int myid, sleeptime;
#ifdef MPI
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);
#else
  myid=ga_nodeid_();
#endif
  sleeptime=(myid+1)/((long) *factor);
#ifdef DEBUG
  printf(" %d sleeping for %d usec %d factor \n", myid, sleeptime, (long) *factor);
#endif
  usleep(sleeptime);
}
