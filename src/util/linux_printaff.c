#define _GNU_SOURCE
#include <sched.h>
#include <stdio.h>
#include <unistd.h>
#ifndef __cpu_set_t_defined
int linux_printaff_(){
  return 0;
}
int linux_setffaff_(){
  return 0;
}
#else
#ifdef MPI
#include <mpi.h>
#else
#include "ga.h"
#include "macdecls.h"
#endif
#define MXCPUS 64
unsigned int i, caff[MXCPUS], numaff=0;
cpu_set_t mycpuid;
int myrank;
pid_t mypid;
int linux_printaff_(){
  mypid=getpid();
#ifdef MPI
  MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
#else
  myrank=GA_Nodeid();
#endif
  CPU_ZERO(&mycpuid);
  if (sched_getaffinity(mypid, sizeof(mycpuid), &mycpuid) < 0) {
    perror("sched_getaffinity");
    return -1;
  }
  for (i = 0; i < MXCPUS; i++){
    if (CPU_ISSET(i, &mycpuid)) {
     caff[numaff]=i;  numaff+=1;
    }
  }
  if(numaff>0) {
    printf("rank %d pid %d bind to %d CPUs:", myrank, (int) mypid, numaff);
    for (i = 0; i < numaff; i++){
      printf(" %i ", caff[i]);
    }
    printf(" \n");
    fflush(stdout);
  }
  return 0;
}
int linux_unsetaff_(){
  mypid=getpid();
#ifdef MPI
  MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
#else
  myrank=GA_Nodeid();
#endif
  CPU_ZERO(&mycpuid);
  for (i = 0; i < MXCPUS; i++){
    CPU_SET(i, &mycpuid);
    }

  if (sched_setaffinity(mypid, sizeof(mycpuid), &mycpuid) < 0) {
    perror("sched_getaffinity");
    return -1;
  }
  return 0;
}
#endif
/* $Id$ */
