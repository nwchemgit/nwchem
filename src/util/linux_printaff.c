#include "typesf2c.h"
#ifndef linux
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
#endif
#define __USE_GNU

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <sched.h>
#ifdef USE_OPENMP
#include <omp.h>
#endif
#include <sys/syscall.h>
#define MXCPUS 1024
unsigned int i, caff[MXCPUS], numaff=0;
cpu_set_t mycpuid;
int myrank, thread;
char hname[128];
Integer linux_printaff_(){
  char affbuf[200 + 7*CPU_SETSIZE];
#ifdef MPI
  MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
#else
  myrank=GA_Nodeid();
#endif
  memset(affbuf,0, sizeof(affbuf));
  memset(hname, 0, sizeof(hname));
  (void)gethostname(hname, sizeof(hname));
#ifdef USE_OPENMP
  #pragma omp parallel private(thread, mycpuid, caff, numaff, affbuf, i)
#endif
  {
    numaff=0;
    CPU_ZERO(&mycpuid);
    (void) sched_getaffinity(0, sizeof(mycpuid), &mycpuid);
#ifdef USE_OPENMP
    thread = omp_get_thread_num();
#else
    thread = 0;
#endif
  for (i = 0; i < MXCPUS; i++){
    if (CPU_ISSET(i, &mycpuid)) {
     caff[numaff]=i;  numaff+=1;
    }
  }
  if(numaff>0) {
    sprintf(affbuf,"rank %d thread %d host %s bind to %d CPUs:", myrank,  thread, hname, numaff);
    for (i = 0; i < numaff; i++){
      sprintf(affbuf + strlen(affbuf)," %i ", caff[i]);
    }
  }
#ifdef USE_OPENMP
#pragma omp barrier
#endif
  puts(affbuf);
  fflush(stdout);
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
