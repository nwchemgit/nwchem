#ifndef linux
#warning __cpu_set_t_defined
Integer linux_printaff_(){
  return (Integer) 0;
}
int linux_setffaff_(){
  return (Integer) 0;
}
#else
#ifdef MPI
#include <mpi.h>
#else
#include "ga.h"
#include "macdecls.h"
#include "typesf2c.h"
#endif
#define __USE_GNU

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <sched.h>
#include <omp.h>
#include <sys/syscall.h>
#define MXCPUS 1024
unsigned int i, caff[MXCPUS], numaff=0;
cpu_set_t mycpuid;
int myrank, thread;
char hname[128];
Integer linux_printaff_(){
#ifdef MPI
  MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
#else
  myrank=GA_Nodeid();
#endif
  memset(hname, 0, sizeof(hname));
  (void)gethostname(hname, sizeof(hname));
  #pragma omp parallel private(thread, mycpuid, caff, numaff)
  {
    numaff=0;
    CPU_ZERO(&mycpuid);
    (void) sched_getaffinity(0, sizeof(mycpuid), &mycpuid);
    thread = omp_get_thread_num();
  for (i = 0; i < MXCPUS; i++){
    if (CPU_ISSET(i, &mycpuid)) {
     caff[numaff]=i;  numaff+=1;
    }
  }
  if(numaff>0) {
    printf("rank %d thread %d host %s bind to %d CPUs:", myrank,  thread, hname, numaff);
    for (i = 0; i < numaff; i++){
      printf(" %i ", caff[i]);
    }
    printf(" \n");
    fflush(stdout);
  }
#pragma omp barrier
  }
  return (Integer) 0;
}
Integer linux_unsetaff_(){
#ifdef MPI
  MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
#else
  myrank=GA_Nodeid();
#endif
  CPU_ZERO(&mycpuid);
  for (i = 0; i < MXCPUS; i++){
    CPU_SET(i, &mycpuid);
    }

  if (sched_setaffinity(0, sizeof(mycpuid), &mycpuid) < 0) {
    perror("sched_getaffinity");
    return (Integer) -1;
  }
  return (Integer) 0;
}
#endif
/* $Id$ */
