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
/*
  this is an example of a C driver calling the
  general eigensystem solver

  Routine geneig_res is currently used to compute
  the residual A Z - B Z D assuming A and B are
  tridiagonal.  If a general A and/or B are used
  then geneig_res must be replaced by residual
  and copies of A and B need to be reset after
  call pdspgv, but before residual

  */

#include <stdio.h>
#include <memory.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>

#include "globalp.c.h"

static Integer countlist();

#ifdef TIMING
#include "timing.h"
TIMINGG test_timing;
#endif

#define ZERO ((DoublePrecision) 0.0e0)
#define max(a,b) ((a) > (b) ? (a) : (b))
#define min(a,b) ((a) < (b) ? (a) : (b))

void main1_()
{
  
  /*
    for  xpress.com
    */
  
  static Integer IONE = (Integer) 1;
  static Integer IZERO = (Integer) 0;
  Integer index;

  Integer n, ii, me, indx, k, i, iii;
  Integer *mapA, *mapB, *mapZ;
  Integer *iscratch;
  DoublePrecision **iptr;
  Integer rsize, ptr_size;
  Integer nprocs, isize;
  Integer info;
  
  DoublePrecision *scratch, *eval, *dptr;
  DoublePrecision *diagA, *subdiagA, *diagB, *subdiagB;
  DoublePrecision *matrixA, *matrixB, *matrixZ;
  DoublePrecision **vecA, **vecB, **vecZ;
  DoublePrecision res, t_com;

#ifdef TIMING
  DoublePrecision time1;
  extern TIMINGG test_timing;
#endif
  
  extern void geneig_res();
  extern void b_ortho();
  extern void tim_com();
  extern void mxend_();
  extern void mxinit_(), mxtime_();
  extern void mxpara_();
  extern DoublePrecision mxclock_(); 
  extern Integer mxnprc_();
  extern Integer mxmynd_();
  extern Integer mxsync_();
  
  extern void memreq_();
  
  extern Integer ci_size_();
  extern DoublePrecision dasum_();

  extern void pdspgv();

  mxinit_();
  me = mxmynd_();
  nprocs = mxnprc_();
  
#ifdef TIMING
  test_timing.choleski = 0.0e0;
  test_timing.inverse  = 0.0e0;
  test_timing.conjug  = 0.0e0;
  test_timing.householder  = 0.0e0;
  test_timing.pstebz  = 0.0e0;
  test_timing.pstein  = 0.0e0;
  test_timing.mxm5x  = 0.0e0;
  test_timing.mxm25  = 0.0e0;
  test_timing.pdspevx  = 0.0e0;
  test_timing.pdspgvx  = 0.0e0;
#endif

  k = 0;
  n = 1020;
  
  diagA    = (DoublePrecision *) malloc( n * sizeof(DoublePrecision));
  subdiagA = (DoublePrecision *) malloc( n * sizeof(DoublePrecision));
  diagB    = (DoublePrecision *) malloc( n * sizeof(DoublePrecision));
  subdiagB = (DoublePrecision *) malloc( n * sizeof(DoublePrecision));

  if( diagA == NULL || subdiagA == NULL || diagB == NULL || subdiagB == NULL ) {
    fprintf(stderr, " me = %d: ERROR not enough memory for diagA or subdiagA, ...\n",
            me );
    exit(-1);
  }
  
  iscratch = (Integer *) malloc ( (4*n + 100) * sizeof(Integer));
  
  if ((mapA = (Integer *) malloc( n * sizeof(Integer))) == NULL ) {
    fprintf(stderr, " me = %d: ERROR not enough memory for mapA %d \n", me, n  );
    exit(-1);
  }
  
  if ((mapB = (Integer *) malloc( n * sizeof(Integer))) == NULL ) {
    fprintf(stderr, "me = %d: ERROR in memory allocation, not enough memory for mapB \n", me);
    exit(-1);
  }
  
  if ((mapZ = (Integer *) malloc( n * sizeof(Integer))) == NULL ) {
    fprintf(stderr, " ERROR in memory allocation, not enough memory for mapZ \n");
    exit(-1);
  }
  
  /*
     set the column mapping of processors
     */
  
  for ( ii = 0;  ii < n; ii++ ) {
    indx = ( ii % nprocs);
    
/*
 *   if ( indx == 3 )
 *     indx = 1;
 *   if ( indx == 2 )
 *     indx = 0;
 */

    
    mapA[ii] = indx;
    mapB[ii] = indx;
  }
  
  for ( ii = 0 ;  ii < n; ii++ ) {
    indx = ( ii % nprocs);
/*
 *   if ( indx == 3 )
 *     indx = 1;
 *   if ( indx == 2 )
 *     indx = 0;
 */
    
    mapZ[ii] = indx;
  }
  
  
  /*
     if ( nprocs > 2 ) {
     mapZ[0] = nprocs-1;
     for ( ii = 1; ii < n; ii++) {
     indx = ( ii  % (nprocs - 1));
     mapZ[ii] = indx;
     }
     }
     else {
     for ( ii = 0; ii < n; ii++) {
     indx = ( ii  % nprocs );
     mapZ[ii] = indx;
     }
     }
     */
  
  /*
     for symmetric matrix with this data distribution
     */

  
  ii = ci_size_( &me, &n, mapA );
  if ( ii > 0 ) {
    if ( (matrixA = (DoublePrecision *) malloc( ii * sizeof(DoublePrecision))) == NULL ) {
      fprintf(stderr, " me %d ERROR in memory allocation, not enough memory for matrixA memory size = %d \n", me, ii);
      exit(-1);
    }
  }
  
  dptr = matrixA;
  for ( indx = 0; indx < ii; indx++ ) {
    *( dptr++ ) = 0.0e0;
  }
  
  ii = countlist ( me, mapA, &n );
  if ( ii > 0 ) {
    if ( ( vecA = ( DoublePrecision ** ) malloc ( ii * sizeof(DoublePrecision *))) == NULL ) {
      fprintf(stderr, "me = %d: ERROR in memory allocation, not enough memory for vecA %d \n", me, ii );
      exit(-1);
    }
  }
  else {
    if ( ( vecA = ( DoublePrecision ** ) malloc ( n * sizeof(DoublePrecision *))) == NULL ) {
      fprintf(stderr, "me = %d: ERROR in memory allocation, not enough memory for vecA %d \n", me, n );
      exit(-1);
    }
  }
  
      
  i = 0;
  dptr = matrixA;
  for ( indx = 0; indx < n; indx++ ) {
    if ( mapA[indx] == me ) {
      vecA[i] = dptr;
      i++;
      dptr += ( n - indx);
    }
  }
  
  i = 0;
  for ( indx = 0; indx < n; indx++ ){

    /*
     * A is symmetric, tri-diagonal.  Set diagA, subdiagA equal
     * to diagonal and subdiagonal parts of matrix.
     * diagA and subdiagA are used to compute residual.
     */

    diagA[indx]    = 1.0/( indx + 1 );
    subdiagA[indx] = -1.0e0;

    if ( mapA[indx] == me ) {
      vecA[i][0] = 1.0/( indx + 1 );
      if ( indx != (n-1))
	vecA[i][1] = -1.0e0;
      i++;
    }
  }
  subdiagA[0] = 0.0e0;

  ii = ci_size_( &me, &n, mapB );
  if ( (matrixB = (DoublePrecision *) malloc( ii * sizeof(DoublePrecision))) == NULL ) {
    fprintf(stderr, " me %d ERROR in memory allocation, not enough memory for matrixB \n", me);
    exit(-1);
  }
  
  zero_out ( ii, matrixB);
  dptr = matrixB;
  for ( indx = 0; indx < ii; indx++ ) {
    *( dptr++ ) = 0.0e0;
  }
  
  ii = countlist ( me, mapB, &n );
  if ( ( vecB = ( DoublePrecision ** ) malloc ( ii * sizeof(DoublePrecision *))) == NULL ) {
    fprintf(stderr, "me = %d: ERROR in memory allocation, not enough memory for vecA \n", me);
    exit(-1);
  }
  
  
  i = 0;
  dptr = matrixB;
  for ( indx = 0; indx < n; indx++ ) {

    /*
     * B is symmetric, tri-diagonal.  Set diagB, subdiagB equal
     * to diagonal and subdiagonal parts of matrix.
     * diagB and subdiagB are used to compute residual.
     */

    diagB[indx]    = 20.0e0;
    subdiagB[indx] = -1.0e0;

    if ( mapB[indx] == me ) {    /* column */
      vecB[i] = dptr;
      vecB[i][0] = 20.0e0;
      if ( indx != ( n-1))
	vecB[i][1]= -1.0e0;
      dptr += ( n-indx);
      i++;
    }
  }
  subdiagB[0] = 0.0e0;
  
  /*
    use the utility routine count_list to determine the number of columns of Z that are stored
    on this processor using the above distribution
    */
  
  ii = countlist ( me, mapZ, &n );
  if ( ( vecZ = ( DoublePrecision ** ) malloc ( ii * sizeof(DoublePrecision *))) == NULL ) {
    fprintf(stderr, "me = %d: ERROR in memory allocation, not enough memory for vecA  allocation = %d \n", me, ii);
    exit(-1);
  }
  
  if ( (matrixZ = (DoublePrecision *) malloc( ii * n * sizeof(DoublePrecision))) == NULL ) {
    fprintf(stderr, "me = %d: ERROR in memory allocation, not enough memory for matrixZ \n", me);
    exit(-1);
  }
  
  dptr = matrixZ;
  
  i = ii*n;
  zero_out( i, matrixZ );
  
  dptr = matrixZ;
  k = 0;
  for ( i = 0; i < ii; i++ ) {
    vecZ[i] = dptr;
    dptr += n;
  }
  
  if ( (eval = (DoublePrecision *) malloc( n * sizeof(DoublePrecision ))) == NULL ) {
    fprintf(stderr, "me = %d: ERROR in memory allocation, not enough memory for eigenvalue space \n", me);
    exit(-1);
  }
  
  index = 0;
  
/*
 * fprintf(stderr, "me = %d: just before memreq \n", me);
*/
  
  rsize = 0;
  isize = 0;
  ptr_size = 0;

/*
  for ( iii = 0; iii < n; iii++ )
    fprintf(stderr, " me = %ld iii = %ld mapA = %ld mapB = %ld mapZ = %ld \n", me, iii, mapA[iii], mapB[iii], mapZ[iii]);
  
*/

  memreq_( &index, &n, mapA, mapB, mapZ, &isize, &rsize, &ptr_size, iscratch );

/*
 * fprintf(stderr, "me = %d: just after memreq isize = %d rsize = %d ptr_size %d \n", me, isize, rsize, ptr_size);
*/
  
  free(iscratch);

  if ( (iscratch = (Integer *) malloc( 2*isize * sizeof(Integer))) == NULL ) {
    fprintf(stderr, " me = %d ERROR in memory allocation, not enough memory for integer scratch space \n", me);
    exit(-1);
  }
  
  rsize = 2 * rsize;
  if ( (scratch = (DoublePrecision *) malloc( (rsize+n+n) * sizeof(DoublePrecision))) == NULL ) {
    fprintf(stderr, " me %d  ERROR in memory allocation, not enough memory for DoublePrecision scratch space \n", me);
    exit(-1);
  }
  
  
  if ( (iptr = (DoublePrecision **) malloc( 2*ptr_size * sizeof(DoublePrecision *))) == NULL ) {
    fprintf(stderr, " me %d ERROR in memory allocation, not enough memory for pointer scratch space \n", me);
    exit(-1);
  }
  
  mxsync_();
  
  if( me == 0 )
    fprintf(stderr, " geneig \n" );
  
#ifdef TIMING
  mxsync_();
  time1 = mxclock_();
#endif

  indx = 1;  
  for ( iii = 0; iii < 1; iii++ ){
    mxtime_( &IZERO, &t_com );
    printf(" in pdspgv node = %d \n", me);
    pdspgv ( &indx, &n, vecA, mapA, vecB, mapB, vecZ, mapZ, eval, iscratch,
	    &isize, iptr, &ptr_size ,scratch, &rsize, &info);
  }
  
  fflush(stdout);
  
#ifdef TIMING
  mxsync_();
  test_timing.pdspgvx = mxclock_() - time1;

  mxtime_( &IONE, &t_com );
  
  ii = 0;
  if ( n < 30 ){
    if ( info == 0 ) {
      for ( k = 0; k < n; k++ ) {
	if ( mapZ[k] == me )  {
	  *scratch = dasum_( &n , vecZ[ii], &IONE );
	  ii++;
	}
      }
    }
  }
  
  if (me == 0 ){
    fprintf(stderr, " n = %d nprocs = %d \n", n, nprocs);
    fprintf(stderr, " pdspgvx = %f \n", test_timing.pdspgvx);
    fprintf(stderr, " pdspevx = %f \n", test_timing.pdspevx);
    fprintf(stderr, " choleski = %f \n", test_timing.choleski);
    fprintf(stderr, " inverse = %f \n", test_timing.inverse);
    fprintf(stderr, " conjug = %f \n", test_timing.conjug);
    fprintf(stderr, " householder = %f \n", test_timing.householder);
    fprintf(stderr, " mxm5x = %f \n", test_timing.mxm5x);
    fprintf(stderr, " mxm25 = %f \n", test_timing.mxm25);
    fprintf(stderr, " pstein = %f \n", test_timing.pstein);
    fprintf(stderr, " pstebz = %f \n", test_timing.pstebz);
  }

  /*  Compute and print commmunication time */

  /*
  tim_com( test_timing.pdspgvx, t_com, iscratch, scratch );
  */

#endif

  geneig_res( &n, diagA, subdiagA, diagB, subdiagB, vecZ, mapZ, eval,
              iscratch, scratch, &res, &info);

  if (me == 0 )
    fprintf(stderr, " A Z - B Z D residual = %g \n", res);

  i = 0;
  for ( indx = 0; indx < n; indx++ ) {
    if ( mapB[indx] == me ) {
      ii = n-indx;
      zero_out( ii, vecB[i] );

      vecB[i][0] = 20.0e0;
      if ( indx != ( n-1))
	vecB[i][1]= -1.0e0;
      i++;
    }
  }

  mxsync_();
  
  
  if ( ! NO_EVEC )
    b_ortho( &n, vecB, mapB, &n, vecZ, mapZ, iptr, iscratch, scratch, &res, &info);
  
  if( me == 0 )
    printf(" (Z' B Z - I)/ulp residual = %g \n", res);
  
  ii = 0;
  if ( n < 30 ){
    if ( info == 0 ) {
      for ( k = 0; k < n; k++ ) {
	if ( mapZ[k] == me )  {
	  *scratch = dasum_( &n , vecZ[ii], &IONE );
	  ii++;
	}
      }
    }
  }
  
  free(iptr);
  free(scratch);
  free(iscratch);
  free(eval);
  free(matrixZ);
  free(vecZ);
  free(vecB);
  free(matrixB);
  free(vecA);
  free(matrixA);
  free(mapZ);
  free(mapB);
  free(mapA);

  return;
  
  /*
     mxpend_();
     */
  
}

static Integer countlist ( me, list, size )
     Integer me, *list, *size;
{
  /*
     
     count the number of instance of "me" in a "list"
    
     */
  
  Integer i, j;
  Integer *ptr;
  
  ptr = list;
  j = 0;
  
  if ( *size <= 0 )
    return(0);
  
  for ( i = 0; i < *size; i++ ) {
    if ( *(ptr++) == me  ) 
      j++;
  }
  return(j);
}

void geneig_res( n, dA, eA, dB, eB, colZ, mapZ, eval, iwork, work, res, info)
     Integer *n, *mapZ, *iwork, *info;
     DoublePrecision *dA, *eA, *dB, *eB, **colZ, *eval, *work, *res;
     
/*
       Given TRIDIAGONAL A and B this routine computes the residual

       res = max_{i = 1 to n} | A z_i - lambda_i * B * z_i | / ( ulp * normA )

       where |.| is the vector infinity-norm and normA is the matrix
       1-norm of A.

       dA[0:n-1] = (input) diagonal of A
       eA[0]     = junk
       eA[1:n-1] = (input) sub-diagonal of A
                        = super-diagonal of A

       dB[0:n-1] = (input) diagonal of B
       eB[0]     = junk
       eB[1:n-1] = (input) sub-diagonal of B
                        = super-diagonal of B
*/

{
  
  Integer           ll, i, j, nvecsZ, me, nprocs, k;
  Integer          *iscrat, *proclist;
  DoublePrecision   t, derror, normA, ulp, t1, t2; 
  DoublePrecision  *ptr, *scrat;
  
  extern DoublePrecision  dlamch_();
  extern Integer          mxmynd_();
  extern Integer          count_list(), reduce_list2();
  extern void             gmax00();

  ll = *n;
  
  *res = 0.e0;

  if ( ll < 1 )
    return;

  me = mxmynd_();
  
  iscrat = iwork;  
  scrat  = work;

  nvecsZ = count_list( me, mapZ, n );
  
  if( nvecsZ == 0 )
    return;
  
  proclist = iscrat;
  nprocs = reduce_list2( *n, mapZ, proclist );

  k      = 0;
  derror =  0.0e0;
  for ( i = 0; i < *n; i++ ) {

    if ( mapZ[i] == me ) {

      ptr = colZ[k];

      if( ll > 1 ) {

        t1 = dA[0] * ptr[0] + eA[1] * ptr[1];
        t2 = dB[0] * ptr[0] + eB[1] * ptr[1];
        t  = fabs( t1 - eval[i] * t2 );

        t1 = dA[ll-1] * ptr[ll-1] + eA[ll-1] * ptr[ll-2];
        t2 = dB[ll-1] * ptr[ll-1] + eB[ll-1] * ptr[ll-2];
        t  = max( t, fabs( t1 - eval[i] * t2 ) );

        for ( j = 1; j < ll-1; j++ ) {
           t1 = eA[j] * ptr[j-1] + dA[j] * ptr[j] + eA[j+1] * ptr[j+1];
           t2 = eB[j] * ptr[j-1] + dB[j] * ptr[j] + eB[j+1] * ptr[j+1];
           t  = max( t, fabs( t1 - eval[i] * t2 ));
        }

      } else {

        t  = fabs( dA[0] * ptr[0] - eval[i] * dB[0] * ptr[0] );

      } 

      derror = max( t, derror);

      k++;
    }
  }
  
  gmax00( (char *) &derror, 1, 5, 16, proclist[0], nprocs, proclist, scrat);
  
      
  if( ll == 1 ) {

    normA = fabs( dA[0] );

  } else {

    normA = fabs( dA[0] ) + fabs( eA[1] );
    normA = max( normA, (fabs( dA[ll-1] ) + fabs( eA[ll-1] )) );
    for ( j = 1; j < ll-1; j++ )
      normA = max( normA, (fabs( dA[j] ) + fabs( eA[j] ) + fabs( eA[j+1] )) );
  }

  if( normA == 0.0e0 )
    normA = 1.0e0;

  ulp = dlamch_( "epsilon") * dlamch_( "base");

  *res = derror / normA / ulp;
  
  return;
}
