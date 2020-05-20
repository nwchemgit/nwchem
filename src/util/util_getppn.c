/* $Id$ */
/* computes the number of processes per node a.k.a ppn 
 it is called by every process only once,
 later calls might not be collective */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <mpi.h>
#include "ga.h"
#include "ga-mpi.h"
#include "typesf2c.h"

#if defined(__bgq__)
#include <process.h>
#include <location.h>
#include <personality.h>
#elif defined(__CRAYXT) || defined(__CRAYXE) || defined(__CRAYXC)
#include <pmi.h>
#endif

static inline int util_mpi_check(int rc, char * name)
{
  if (rc != MPI_SUCCESS) {
    fprintf(stdout,"util_getppn: %s failed\n",name);
    return 1;
  }
  return 0;
}

static short int ppn_initialized=0;
static int ppn=0;
void FATR util_getppn_(Integer *ppn_out){

#if defined(__bgq__)
  *ppn_out = (Integer) Kernel_ProcessCount();
  return;
  if(0) {
#elif MPI_VERSION >= 3

    int err;
    MPI_Comm comm_node;

  if(ppn_initialized) {
    *ppn_out = (Integer) ppn;
    
  }else{
    MPI_Comm ga_comm=GA_MPI_Comm_pgroup_default();
    err = MPI_Comm_split_type(ga_comm, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &comm_node);
    if (util_mpi_check(err,"MPI_Comm_split_type")) goto errlab;

    err = MPI_Comm_size(comm_node, &ppn);
    if (util_mpi_check(err,"MPI_Comm_size")) goto errlab;

    err = MPI_Comm_free(&comm_node);
    if (util_mpi_check(err,"MPI_Comm_free")) goto errlab;

    ppn_initialized=1;
    *ppn_out = (Integer) ppn;
    return;
#else // no MPI-3 or machine-specific optimized implementation
    /* A space-efficient implementation is described in pseudo-code here:
     * http://lists.mcs.anl.gov/pipermail/mpich-discuss/2012-January/011662.html
     * A slightly more efficient implementation than below may be:
     * https://github.com/jeffhammond/HPCInfo/blob/master/mpi/advanced/symmetric-heap.c */
  const int mxlen = 255;
  char myhostname[mxlen];
  char* recvbuf;
  int i, num_procs, me,  err, modppn;
  MPI_Comm ga_comm=GA_MPI_Comm_pgroup_default();
  
  if(ppn_initialized) {
    *ppn_out = (Integer) ppn;
    
  }else{
    num_procs = GA_Nnodes();
    me = GA_Nodeid();
    
      recvbuf=(char*)malloc(num_procs*(mxlen+1)*(sizeof(char)));
      
      err=gethostname(myhostname, sizeof(myhostname) );
      if (err != 0) {
	fprintf(stdout,"util_getppn: gethostname failed\n");
	ppn=0;
	goto errlab;
      }
      
      
      err=MPI_Allgather(myhostname, mxlen, MPI_CHAR, recvbuf, mxlen, MPI_CHAR, ga_comm);
      if (err != MPI_SUCCESS) {
	fprintf(stdout,"util_getppn: MPI_Allgather failed\n");
	ppn=0;
	goto errlab;
      }
      
      
      for (i=0; i< num_procs; i++){
	if(strcmp(myhostname,&recvbuf[mxlen*i])==0) ppn++;
      }
      
      /*	  free malloc'ed memory */
      free(recvbuf);
      
      
    /* broadcast ppn to everybody */
    err= MPI_Bcast(&ppn, 1, MPI_INT, 0, ga_comm);
    if (err != MPI_SUCCESS) {
      fprintf(stdout,"util_getppn: MPI_Bcast failed\n");
      goto errlab;
    }
    

    /* check that computed ppn is a submultiple of num procs */
    
    modppn = num_procs%ppn;
    if (modppn ==0){
      ppn_initialized=1;
      *ppn_out = (Integer) ppn;
      return;
    }else{
      printf(" ERROR: numprocs %d  ppn %d  mod %d\n", num_procs, ppn,  modppn);
      goto errlab;
    }
#endif
  errlab:
    GA_Error(" util_getppn failure", 0);
    return;
  }
} 

/* C binding for util_ppn */
int util_cgetppn(){
  Integer* ppn_out = malloc(sizeof(Integer));
  int ppn;
  util_getppn_(ppn_out);
  ppn = (int ) *ppn_out;
  free(ppn_out);
  fflush(stdout);
  return ppn;
}



int util_my_smp_index(){
  int ppn= util_cgetppn();
  return GA_Nodeid()%ppn;
}

int util_my_smp_master(){
  int ppn= util_cgetppn();
  return (GA_Nodeid()/ppn)*ppn;
}

