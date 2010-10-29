/*
 $Id$
 *======================================================================
 *
 * DISCLAIMER
 *
 * This material was prepared as an account of work sponsored by an
 * agency of the United States Government.  Neither the United States
 * Government nor the United States Department of Energy, nor Battelle,
 * nor any of their employees, MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR
 * ASSUMES ANY LEGAL LIABILITY OR RESPONSIBILITY FOR THE ACCURACY,
 * COMPLETENESS, OR USEFULNESS OF ANY INFORMATION, APPARATUS, PRODUCT,
 * SOFTWARE, OR PROCESS DISCLOSED, OR REPRESENTS THAT ITS USE WOULD NOT
 * INFRINGE PRIVATELY OWNED RIGHTS.
 *
 * ACKNOWLEDGMENT
 *
 * This software and its documentation were produced with Government
 * support under Contract Number DE-AC06-76RLO-1830 awarded by the United
 * States Department of Energy.  The Government retains a paid-up
 * non-exclusive, irrevocable worldwide license to reproduce, prepare
 * derivative works, perform publicly and display publicly by or for the
 * Government, including the right to distribute to other Government
 * contractors.
 *
 *======================================================================
 *
 *  -- PEIGS  routine (version 2.1) --
 *     Pacific Northwest Laboratory
 *     July 28, 1995
 *
 *======================================================================
 */
#include <stdio.h>

#include "globalp.c.h"

void ortho2( n, m, colZ, mapZ, ibuffptr, iwork, work, ort, info)
     Integer *n, *m, *mapZ, *iwork, *info;
     DoublePrecision **colZ, **ibuffptr, *work,  *ort;
     
     /*
       
   This subroutine computes the infinity-norm measure:
   
   ort = max_i | (Z^t.Z)_i - I_i | / ULP,
   
   for the standard symmetric eigensystem problem where
   
   Z is N-by-M
   I is M-by-M.
   Z_i is an eigenvector,
   (Z^t.Z)_i is the i-th column of Z^t.Z
   I_i is the i-th column of the m-by-m identity matrix
   ULP = (machine precision) * (machine base)
   |.| is the infinity-norm.
   
   res is reasonable if it is of order about 50 or less.

   
   MUST have M <= N.  If M > N then program exits.
   This is not a limitation of this subroutine as M > N
   implies the columns of Z are linearly dependent, which
   implies "ort" will always be large in this case.

   Arguments
   ---------
   In the following:
   
   INTEGER          = "pointer to Integer"
     DOUBLE PRECISION = "pointer to DoublePrecision"
     
     me     = this processor's id (= mxmynd_())
     nprocs = number of allocated processors ( = mxnprc_())
     nvecsZ = number of entries in mapZ equal to me
     (= count_list( me, mapZ, m ))
     nvecsZ_max = maximum number of entries in mapZ equal to i,
     i = 0 to nprocs-1
     (= max over i of count_list( i, mapZ, m ))
     sDP    = sizeof( DoublePrecision )
     
     
     n....... (input) INTEGER
     Number of rows in Z
     
     m....... (input) INTEGER
     number of columns in Z (i.e., # of 
            eigenvalues/eigenvectors).
            Must have m <= n.
	    
	    colZ ... (input) array of pointers to DoublePrecision,
	    length (nvecsZ)
            The part of matrix Z owned by this processer, stored
            in packed format, i.e., colZ[i] points to the start
            of the i-th column of matrix Z owned by this
            processor, i = 0 to nvecsZ-1.
                
   mapZ ... (input) INTEGER array, length (m)
            The i-th column of matrix Z is owned by processor
            mapZ[i], i = 0 to m-1.

   ibuffptr (workspace) array of pointers to DoublePrecision,
                        length(nvecsZ)

   iwork... (workspace) INTEGER array, length( 7 * m )

   work.... (workspace) DOUBLE PRECISION array,
                        length( nvecsZ * n + maximum( d_mxm25,
                                                      mxlbuf_() / sDP + 1 )
                        where d_mxm25 = maximum ( 2*m, (nvecsZ + 2*nvecsZ_max)*n )
       
   ort..... (output) INTEGER
            the residual described above.

   info.... (output) INTEGER
            = 0, not currently used
       
 */
{
  
  static Integer ONE = 1;
  
  Integer ll, i, j, *iscrat, *mapvecZ;
  Integer nvecsZ;
  Integer me;
  Integer nprocs, *proclist;
  
  DoublePrecision t, ulp;
  DoublePrecision *ptr, *scrat;
  DoublePrecision **vecZ1; /* copies of the vecZ matrix */
  
  /*
    blas call
    */
  
  extern void dcopy_();
  extern void mxm25();
  
  extern Integer mxmynd_();

  extern Integer fil_mapvec_();
  extern Integer reduce_list2();
  extern void mat_max();
  
  /*
    usual story about error handling
    */
  
  ll = *n;

  *ort = 0.e0;
  *info = 0;
  
  if ( ll < 1 )
    return;  /* error */

  me = mxmynd_();
  
  if( *n < *m ) {
      fprintf(stderr,
              "Error in routine ortho m (=%d) < n (=%d). me = %d \n",
              *m, *n, me);
      exit(-1);
  }
  
  iscrat = iwork;
  mapvecZ = iscrat;
  nvecsZ = fil_mapvec_( &me, m, mapZ, mapvecZ );
  iscrat += nvecsZ;

  proclist = iscrat;
  nprocs = reduce_list2( *m, mapZ, proclist);
  iscrat += nprocs;
  
  if ( nvecsZ == 0 )
    return;
  
  scrat = work;
  
  vecZ1 = ibuffptr;
  
  /*
    copy over the matrix to a buffer area
    */
  
  ptr = scrat;
  for ( i = 0; i < nvecsZ; i++ ) {
    printf(" i = %d \n", i);
    vecZ1[i] = ptr;
    dcopy_( &ll, colZ[i], &ONE, ptr, &ONE );
    ptr += ll;
  }
  
  scrat = ptr;
  
  
  /*
    Z1 is a copy of Z
    
    Z1 <- Z^t . Z1;
  */
  
  mxm25 ( m, &ll, colZ, mapZ, m, vecZ1, mapZ, vecZ1, iscrat, scrat);
  
  for ( i = 0; i < nvecsZ; i++ ) {
    j = mapvecZ[i];
    vecZ1[i][j] -= 1.0e0;
  }
  
  mat_max ( m, m, vecZ1, mapZ, ort, iscrat, scrat);
  
  return;
}
