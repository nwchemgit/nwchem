/* $Id$ */
/* routine to avoid 32-bit integer overflow present both in GA and MPI collectives*/
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include "ga.h"
#include "macdecls.h"
#include "typesf2c.h"

void FATR
util_mygabcast_(Integer *g_a, Integer *m, Integer *n, DoublePrecision *a, Integer *ld) {
  int i;
  int ierr, len,  resultlen;
  long  len8,  istart, end8;
  int nsteps;
  int64_t lo[2], hi[2], *ld64=NULL;
  static int bigint = 2147482574;
  char err_buffer[MPI_MAX_ERROR_STRING];
  lo[0]=0;
  lo[1]=0;
  /* swap column & rows when going from fortran to c */
  hi[1]=*m - 1;
  hi[0]=*n - 1;

  if(GA_Nodeid() == 0)  NGA_Get64(*g_a, lo, hi, a, (int64_t *) ld);
  GA_Sync();

  len8= (*m) * (*n);

  nsteps = (int) ceil(((double)len8)/((double)bigint));

  istart=1;

  for (i=0; i < nsteps; i++){
  
    len=bigint;
    end8=istart+len-1;
    if (istart+len-1 > len8) {
      len=((long)(len8 - istart)) + 1;
      end8=istart+len-1;
    }

#ifdef DEBUG   
    if(GA_Nodeid() == 0) printf(" bcast: is %11ld len %11d  end8 %11ld step i %2d of %4d lentot %8ld\n", istart, len, istart+(long)(len-1), i+1, nsteps, len8);
#endif

    ierr= MPI_Bcast(a+istart, len, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD);
    
    if (ierr != MPI_SUCCESS) {
      fprintf(stderr,"util_mygabcast: MPI_Bcast failed\n");
      MPI_Error_string(ierr,err_buffer,&resultlen);
      fprintf(stderr,"%s\n", err_buffer);
      GA_Error("util_mygabcast error", 0L);
    }
    
    istart+=len;
  }
  MPI_Barrier(MPI_COMM_WORLD);


}
