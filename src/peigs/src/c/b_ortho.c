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
#include <stdlib.h>

#define MSG_START 25000
#define ZERO ((DoublePrecision) 0.0e0)

#include "globalp.c.h"

void b_ortho ( n, colB, mapB, m, colZ, mapZ, ibuffptr, iwork, work, ort, info)
     Integer *n, *mapB, *m, *mapZ, *iwork, *info;
     DoublePrecision **colB, **colZ, **ibuffptr, *work,  *ort;
     
/*
   This subroutine computes the  infinity-norm measure:
       
   ort = max_i | (Z^t.B.Z)_i - I_i | / ULP,

   for the generalized symmetric eigensystem problem where

     Z is N-by-M and distributed by columns
     B is N-by-N symmetric and in packed storage
     I is M-by-M.
     Z_i is an eigenvector,
     (Z^t.B.Z)_i is the i-th column of Z^t.B.Z
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
                  (= count_list( me, mapZ, n ))
     nvecsZ_max = maximum number of entries in mapZ equal to i,
                  i = 0 to nprocs-1
                  (= max over i of count_list( i, mapZ, n ))
     sDP    = sizeof( DoublePrecision )

       
   n....... (input) INTEGER
            Number of rows and columns in B, and
            the number of rows in Z

   colB ... (input) array of pointers to DoublePrecision,
                    length (nvecsB)
            The part of matrix B owned by this processer stored
            in packed format, i.e., colB[i] points to the diagonal
            element of the i-th column (or equivalently row) of B
            owned by this processor, i = 0 to nvecsB-1.
                
   mapB ... (input) INTEGER array, length (n)
            The i-th column (or equivalently row) of B is 
            owned by processor mapB[i], i = 0 to n-1.

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
                        length( 3*nvecsZ + 1)

   iwork... (workspace) INTEGER array, 
                        length( n+m + MAX{nprocs, 3*n+2*m+nvecsB+nvecsZ }

   work.... (workspace) DOUBLE PRECISION array,
                        length( nvecsZ * (n+m) + maximum( d_mxm25, d_mxm88,
                                                          mxlbuf_() / sDP + 1 )
                        where d_mxm25 = maximum ( n+m, (nvecsZ + 2*nvecsZ_max)*n )
                              d_mxm88 = maximum ( (nvecsZ+1)*n, mxlbuf_()/sDP + 1 )
       
   ort..... (output) INTEGER
            the residual described above.

   info.... (output) INTEGER
            = 0, not currently used
       
 */
{
  
  static Integer IONE = 1;
  
  Integer ll, i, j, *iscrat, *mapvecB, *mapvecZ;
  Integer nvecsB, nvecsZ;
  Integer me, nprocs;
  
  DoublePrecision ulp;
  DoublePrecision *ptr, *scrat;
  DoublePrecision **vecZ1, **vecZ2; /* copies of the vecZ matrix */
  
  /*
    blas call
    */
  
  extern DoublePrecision dnrm2_();
  extern void mxm88();
  extern void mxm25();
  extern DoublePrecision ddot_();
  extern void dcopy_();
  extern void daxpy_();
  
  extern Integer mxwrit_(), mxread_();
  extern Integer menode_(), mxmynd_();
  extern void    mdiff1_(), bbcast00();

  extern void chol_pipe_bcast();
  extern Integer fil_mapvec_();
  extern Integer reduce_list2();
  extern Integer count_list();
  extern void fil_dbl_lst();
  extern void mat_max();
  
  /*
    usual story about error handling
    */
  
  ll = *n;
  
  *info = 0;
  *ort = 0.e0;

  if ( ll < 1 )
    return;  /* error */
  
  /*
    should perform global sync to check for errors
    */
  
  me = mxmynd_();
  
  if( *n < *m ) {
      fprintf(stderr,
              "Error in routine b_ortho. m (=%d) < n (=%d). me = %d \n",
              *m, *n, me);
      exit(-1);
  }
  
  iscrat = iwork;

  mapvecB = iscrat;
  nvecsB = fil_mapvec_( &me, &ll, mapB, mapvecB);
  iscrat += nvecsB;

  mapvecZ = iscrat;
  nvecsZ = fil_mapvec_( &me, m, mapZ, mapvecZ );
  iscrat += nvecsZ;
  
  if ( ( nvecsB == 0 ) && ( nvecsZ == 0 ))
    return;
  
  scrat = work;
  
  vecZ1 = ibuffptr;
  
  /*
    copy over the matrix to a buffer area
    */
  
  ptr = scrat;
  for ( i = 0; i < nvecsZ; i++ ) {
    vecZ1[i] = ptr;
    dcopy_( &ll, colZ[i], &IONE, ptr, &IONE );
    ptr += ll;
  }
  
  vecZ2 = vecZ1 + nvecsZ;
  for ( i = 0; i < nvecsZ; i++ ) {
    vecZ2[i] = ptr;
    fil_dbl_lst( *m, ptr, 0.0);
    ptr += *m;
  }
  
  scrat = ptr;
  
  /*
    Z1 is a copy of Z
    Z1 <- B.Z1
    */
  
  mxm88( &ll, colB, mapB,  m, vecZ1, mapZ, iscrat, scrat, vecZ1 + 2*nvecsZ + 1);
  
  /*
    Z2 <- Z^t . Z1; Z2 is now Z^t. B. Z
    */
  
  mxm25( m, &ll, colZ, mapZ, m, vecZ1, mapZ, vecZ2, iscrat, scrat);
  
  for ( i = 0; i < nvecsZ; i++ ) {
    j = mapvecZ[i];
    vecZ2[i][j] -= 1.0e0;
  }
  
  *ort = 0.0e0;
  mat_max ( m, m, vecZ2, mapZ, ort, iscrat, scrat);
  
  ulp = DLAMCHE * DLAMCHB;
  
  *ort = *ort / ulp;
  mdiff1_( n, mapB, m, mapZ, iscrat, &nprocs ); 
  
  if( nprocs > 0 ){
    iscrat[nprocs] = mapZ[0];
    nprocs++;
    bbcast00( (char *) ort, sizeof(DoublePrecision), 1, mapZ[0], nprocs, iscrat );
  }
  
  return;
}

