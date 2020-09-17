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
 * Test code for the Wilkinson tridiagonal test matrices
 */

#include <stdio.h>
#include <memory.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>

#include "globalp.c.h"

#ifdef TIMING
#include "timing.h"
TIMINGG test_timing;
#endif


void main1_()
{
  
  /*
    for  xpress.com
    */
  
static Integer IZERO = (Integer) 0;
  Integer index;
  
  Integer n, ii, me, indx, k, i, neleZ, neleA;
  Integer *mapA, mapB[1], *mapZ;
  Integer *iscratch;
  DoublePrecision **iptr;
  DoublePrecision *dd, *ee;
  Integer rsize, ptr_size;
  Integer nprocs, isize;
  Integer info, m;
  
  DoublePrecision *scratch, *eval, *dptr;
  DoublePrecision *matrixA, *matrixZ;
  DoublePrecision **vecA, **vecZ;
  DoublePrecision res, t_com;

#ifdef TIMING
  static Integer IONE = (Integer) 1;
  DoublePrecision time1, time2, timex;
  extern TIMINGG test_timing;
#endif
  
  Integer countlist();
  
  extern void tim_com();
  extern void mxend_();
  extern void mxinit_(), mxtime_();
  extern DoublePrecision mxclock_(); 
  extern void mxpara_();
  extern Integer mxnprc_();
  extern Integer mxmynd_();
  
  extern void memreq_();
  
  extern Integer ci_size_();
  extern DoublePrecision dasum_();

  extern void pdspev();
  extern void tresid(), ortho();

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

  /*
  while (1) {
  */
  {
    m = 250;

  n = 2 * m + 1;

  
  nprocs = mxnprc_();
  printf(" n = %d nprocs = %d \n", n, nprocs);
  
  if ((dd = (DoublePrecision *) malloc( n * sizeof(DoublePrecision))) == NULL ) {
    fprintf(stderr, " me = %d: ERROR in memory allocation, not enough memory for dd %d \n", me, n  );
    exit(-1);
  }

  if ((ee = (DoublePrecision *) malloc( n * sizeof(DoublePrecision))) == NULL ) {
    fprintf(stderr, " me = %d: ERROR in memory allocation, not enough memory for ee %d \n", me, n  );
    exit(-1);
  }

  if ((mapA = (Integer *) malloc( n * sizeof(Integer))) == NULL ) {
    fprintf(stderr, " me = %d: ERROR in memory allocation, not enough memory for mapA %d \n", me, n  );
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
    mapA[ii] = indx;
  }
  
  for ( ii = 0 ;  ii < n; ii++ ) {
    indx = ( ii % nprocs);
    mapZ[ii] = indx;
  }
  
  
  /*
     for symmetric matrix with this data distribution
     */

  
  ii = ci_size_( &me, &n, mapA );
  neleA = ii;
  if ( ii > 0 ) {
    if ( (matrixA = (DoublePrecision *) malloc( ii * sizeof(DoublePrecision))) == NULL ) {
      fprintf(stderr, " me %d ERROR in memory allocation, not enough memory for matrixA memory size = %d \n", me, ii);
      exit(-1);
    }
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
  
  /*
    wilkinson's matrix
    */
  
  ee[0] = 0.e0;
  for ( indx = 1; indx < n; indx++)
    ee[indx] = 1.0e0;
  
  /*
    ee[indx] = 1.e0;
  */
  
  i = 0;
  for ( indx = 0; indx < m; indx++){
    dd[indx] = (DoublePrecision) ( m-indx + pow(10., (DoublePrecision) indx));
    if ( mapA[indx] == me ){
      vecA[i][0] = (DoublePrecision) ( m-indx + pow(10., (DoublePrecision) indx));
      vecA[i][1] = 1.;
      i++;
    }
  }
  dd[m] = 0.e0;
  if ( mapA[m] == me ) {
    vecA[i][0] = 0.;
    vecA[i][1] = 1.;
    i++;
  }
  for ( indx = m+1; indx < n; indx++){
    dd[indx] = (DoublePrecision) ( m-indx - pow(10., (DoublePrecision) indx));
    if ( mapA[indx] == me ){
      vecA[i][0] = (DoublePrecision) ( m-indx - pow(10., (DoublePrecision) indx));
      if ( indx != n-1)
	vecA[i][1] = 1;
      i++;
    }
  }

  /*
    use the utility routine count_list to determine the number of columns of Z that are stored
    on this processor using the cve distribution
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
  
  neleZ = ii*n;
  
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
  
  index = 1;
  /*
  fprintf(stderr, "me = %d: just before memreq \n", me);
  */
  rsize = 0;
  isize = 0;
  ptr_size = 0;
  iscratch = ( Integer *) malloc( 6*n*sizeof(Integer));
  memreq_( &index, &n, mapA, mapB, mapZ, &isize, &rsize, &ptr_size, iscratch );
  /*
  fprintf(stderr, "me = %d: just after memreq isize = %d rsize = %d ptr_size %d \n", me, isize, rsize, ptr_size);
  */
  
  free(iscratch);
  
  if ( (iscratch = (Integer *) malloc( 4* isize * sizeof(Integer))) == NULL ) {
    fprintf(stderr, " me = %d ERROR in memory allocation, not enough memory for integer scratch space \n", me);
    exit(-1);
  }
  
  if ( (scratch = (DoublePrecision *) malloc( 4*rsize * sizeof(DoublePrecision))) == NULL ) {
    fprintf(stderr, " me %d  ERROR in memory allocation, not enough memory for DoublePrecision scratch space \n", me);
    exit(-1);
  }
  
  
  if ( (iptr = (DoublePrecision **) malloc( 4*ptr_size * sizeof(DoublePrecision *))) == NULL ) {
    fprintf(stderr, " me %d ERROR in memory allocation, not enough memory for pointer scratch space \n", me);
    exit(-1);
  }
  
  if( me == 0 )
    fprintf(stderr, " Wilkinson \n" );
  
  for ( ii = 0; ii < 1; ii++ ) {
  
     /* set data modified by pdspevx */
  
     zero_out( neleZ, matrixZ );
     zero_out( neleA, matrixA );
  
     for ( k = 0;  k < n; k++ ) {
       indx = ( k % nprocs);
       mapZ[k] = indx;
     }
  
     k = 0;
     for ( indx = 0; indx < n; indx++ ){
       if ( me == mapA[indx] ) {
         vecA[k][0] = dd[indx];
         if ( indx != n-1 )
   	   vecA[k][1] = ee[indx+1];
         k++;
       }
     }
     
#ifdef TIMING
     time1 = mxclock_();

     mxsync_();

     time1 = mxclock_();
#endif

     mxtime_( &IZERO, &t_com );

     pdspev( &n, vecA, mapA, vecZ, mapZ, eval, iscratch,
	    &isize, iptr, &ptr_size ,scratch, &rsize, &info);

     if ( me == 0 )
       for ( k = 0; k < n; k++ )
	 printf(" driver wilk k = %d eval %f \n", k, eval[k]);
     
     mxsync_();
#ifdef TIMING
     timex = mxclock_();

     mxtime_( &IONE, &t_com );
  
     if( ii == 0 )
       time2 = timex - time1;
#endif
     
     if (!NO_EVEC){
       
       tresid( &n, &n, dd, ee, vecZ, mapZ, eval, iscratch, scratch, &res, &info);
       
       if( me == 0 )
	 fprintf(stderr, " iteration # %d : A Z - Z D residual = %g \n", ii, res);
       
       ortho( &n, &n, vecZ, mapZ, iptr, iscratch, scratch, &res, &info);
       
       if( me == 0 )
	 fprintf(stderr, " iteration # %d : Z' Z - I residual = %g \n", ii, res);
       
     }
  }
     
#ifdef TIMING
  
  test_timing.pdspevx = timex - time1;

  if (!NO_EVEC){
    ii = 0;
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
    fprintf(stderr, " time1   = %f \n", time2);
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

  tim_com( test_timing.pdspevx, t_com, iscratch, scratch );

#endif


  free(iptr);
  free(scratch);
  free(iscratch);
  free(eval);
  free(matrixZ);
  free(vecZ);
  free(vecA);
  free(matrixA);
  free(mapZ);
  free(mapA);
}
  return;
  
  
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

