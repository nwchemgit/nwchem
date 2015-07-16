/* $Id$ */
/* routine to avoid 32-bit integer overflow present both in GA and MPI collectives*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include "ga.h"
#include "macdecls.h"
#include "typesf2c.h"
extern MPI_Comm GA_MPI_Comm_pgroup_default();

void FATR
util_mygabcast2_(Integer *g_a, Integer *mlo, Integer *mhi, Integer *nlo, Integer *nhi, DoublePrecision *a, Integer *ld) {
  int i;
  int ierr, len,  resultlen;
  int nsteps;
  int64_t lo[2], hi[2];
  char err_buffer[MPI_MAX_ERROR_STRING];
  long bigint8 = ((long)pow(2.,31.)-1024)/sizeof(DoublePrecision);
  int bigint = (int) bigint8;

  MPI_Comm ga_comm=GA_MPI_Comm_pgroup_default();

#ifdef DEBUG   
    if(GA_Nodeid() == 0) printf(" bcast: bigint8 %11ld bigint %11d\n", bigint8, bigint);
#endif
  lo[1]=*mlo - 1;
  lo[0]=*nlo - 1;
  /* swap column & rows when going from fortran to c */
  hi[1]=*mhi - 1;
  hi[0]=*nhi - 1;

  long len8 = (hi[1] - lo[1] +1) * (hi[0] - lo[0] +1);

  nsteps = (int) ceil(((double)len8)/((double)bigint));

  long istart=0;

  if(GA_Nodeid() == 0) NGA_Get64(*g_a, lo, hi, a, (int64_t *) ld);
  GA_Sync();

  for (i=0; i < nsteps; i++){
  
    len=bigint;
    if (istart+len-1 > len8) len=((long)(len8 - istart)) + 1;

#ifdef DEBUG   
    if(GA_Nodeid() == 0) printf(" bcast: is %11ld len %11d  step i %2d of %4d lentot %8ld\n", istart, len, istart+(long)(len-1), i+1, nsteps, len8);
#endif

    ierr= MPI_Bcast(a+istart, len, MPI_DOUBLE_PRECISION, 0, ga_comm);
    
    if (ierr != MPI_SUCCESS) {
      fprintf(stdout,"util_mygabcast: MPI_Bcast failed step %2d len %11d\n", i, len);
      MPI_Error_string(ierr,err_buffer,&resultlen);
      fprintf(stdout,"%s\n", err_buffer);
      GA_Error("util_mygabcast error", 0L);
    }
    
    istart+=len;
  }
  MPI_Barrier(ga_comm);


}
#if 1
void FATR
util_mygabcast_(Integer *g_a, Integer *m, Integer *n, DoublePrecision *a, Integer *ld) {
  Integer* mlo = malloc(sizeof(Integer));
  Integer* mhi = malloc(sizeof(Integer));
  Integer* nlo = malloc(sizeof(Integer));
  Integer* nhi = malloc(sizeof(Integer));
  *mlo = 1;
  *mhi = *m;
  *nlo = 1;
  *nhi = *n;
  util_mygabcast2_(g_a, mlo,  mhi,  nlo, nhi, a, ld);
  free(mlo);
  free(mhi);
  free(nlo);
  free(nhi);
 }
#endif
