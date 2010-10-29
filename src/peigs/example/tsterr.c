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
 * Routine to test error returns of drivers.
 *
 * Should delete STOP statement from xstop.f in 
 * ..../src/f77/xstop.f before linking this code.
 */

#include <stdio.h>
#include <memory.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>

#include "globalp.c.h"


#define ZERO ((DoublePrecision) 0.0e0)
#define max(a,b) ((a) > (b) ? (a) : (b))
#define min(a,b) ((a) < (b) ? (a) : (b))

void main1_()
{
  
  Integer index;

  Integer nbad[10], n, ii, me, indx, i;
  Integer *mapA, *mapB, *mapZ;
  Integer *iscratch;
  DoublePrecision **iptr;
  Integer rsize, ptr_size;
  Integer nprocs, isize;
  
  DoublePrecision *scratch, *eval, *dptr;
  DoublePrecision *matrixA, *matrixB, *matrixZ;
  DoublePrecision **vecA, **vecB, **vecZ;

  extern void mxinit_();
  extern Integer mxsync_(), mxnprc_();
  extern Integer mxmynd_();
  
  extern void memreq_();
  
  extern Integer tst_pdspgv(), tst_pdspgvx(), tst_pdspevx(), tst_pdspev();

  mxinit_();
  me = mxmynd_();
  nprocs = mxnprc_();
  
  n = 50;

  n++;
  mapA = (Integer *) malloc( n * sizeof(Integer));
  mapB = (Integer *) malloc( n * sizeof(Integer));
  mapZ = (Integer *) malloc( n * sizeof(Integer));

  if( mapA == NULL || mapB == NULL || mapZ == NULL ) {
    fprintf(stderr, " me = %d: ERROR not enough memory for maps \n", me );
    exit(-1);
  }
  
  /*
   *  set the column mapping of processors
   */
  
  for ( ii = 0;  ii < n; ii++ ) {
    indx = ( ii % nprocs);
    
    mapA[ii] = indx;
    mapB[ii] = indx;
    mapZ[ii] = indx;
  }
  mapA[n-1] = mapA[n-2];
  mapB[n-1] = mapB[n-2];
  mapZ[n-1] = mapZ[n-2];
  
  ii = n * (n+1) / 2;
  matrixA = (DoublePrecision *) malloc( ii * sizeof(DoublePrecision));
  if( matrixA == NULL ) {
    fprintf(stderr, " me = %d: ERROR not enough memory for matrixA \n", me );
    exit(-1);
  }
  
  dptr = matrixA;
  for ( indx = 0; indx < ii; indx++ ) {
    *( dptr++ ) = 0.0e0;
  }
  
  ii = n;
  vecA = ( DoublePrecision ** ) malloc ( ii * sizeof(DoublePrecision *));
  if ( vecA == NULL ) {
    fprintf(stderr, "me = %d: ERROR not enough memory for vecA %d \n", me, ii );
    exit(-1);
  }
      
  i = 0;
  dptr = matrixA;
  for ( indx = 0; indx < n; indx++ ) {
      vecA[i] = dptr;
      i++;
      dptr += ( n - indx);
  }
  
  i = 0;
  for ( indx = 0; indx < n; indx++ ){
      vecA[i][0] = 1.0/( indx + 1 );
      if ( indx != (n-1))
        vecA[i][1] = -1.0e0;
      i++;
  }

  ii = n*(n+1)/2;
  if ( (matrixB = (DoublePrecision *) malloc( ii * sizeof(DoublePrecision))) == NULL ) {
    fprintf(stderr, " me %d ERROR in memory allocation, not enough memory for matrixB \n", me);
    exit(-1);
  }
  
  zero_out ( ii, matrixB);
  dptr = matrixB;
  for ( indx = 0; indx < ii; indx++ ) {
    *( dptr++ ) = 0.0e0;
  }
  
  ii = n;
  if ( ( vecB = ( DoublePrecision ** ) malloc ( ii * sizeof(DoublePrecision *))) == NULL ) {
    fprintf(stderr, "me = %d: ERROR in memory allocation, not enough memory for vecA \n", me);
    exit(-1);
  }
  
  
  i = 0;
  dptr = matrixB;
  for ( indx = 0; indx < n; indx++ ) {
      vecB[i] = dptr;
      vecB[i][0] = 20.0e0;
      if ( indx != ( n-1))
        vecB[i][1]= -1.0e0;
      dptr += ( n-indx);
      i++;
  }
  
  ii = n;
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
  for ( i = 0; i < ii; i++ ) {
    vecZ[i] = dptr;
    dptr += n;
  }
  
  if ( (eval = (DoublePrecision *) malloc( n * sizeof(DoublePrecision ))) == NULL ) {
    fprintf(stderr, "me = %d: ERROR in memory allocation, not enough memory for eigenvalue space \n", me);
    exit(-1);
  }
  
  index = 0;
  
  rsize = 0;
  isize = 0;
  ptr_size = 0;
  iscratch = (Integer *) malloc ( (4*n + 100) * sizeof(Integer));
  memreq_( &index, &n, mapA, mapB, mapZ, &isize, &rsize, &ptr_size, iscratch );
  free(iscratch);

  if ( (iscratch = (Integer *) malloc( 2*isize * sizeof(Integer))) == NULL ) {
    fprintf(stderr, " me = %d ERROR in memory allocation, not enough memory for integer scratch space \n", me);
    exit(-1);
  }
  
  if ( (scratch = (DoublePrecision *) malloc( 2*rsize * sizeof(DoublePrecision))) == NULL ) {
    fprintf(stderr, " me %d  ERROR in memory allocation, not enough memory for DoublePrecision scratch space \n", me);
    exit(-1);
  }
  
  if ( (iptr = (DoublePrecision **) malloc( 2*ptr_size * sizeof(DoublePrecision *))) == NULL ) {
    fprintf(stderr, " me %d ERROR in memory allocation, not enough memory for pointer scratch space \n", me);
    exit(-1);
  }
  
 n--;

 if( me < n ) {
   nbad[0] = tst_pdspgv( n, vecA, mapA, vecB, mapB, vecZ, mapZ, eval,
                          iscratch, isize, iptr, ptr_size ,scratch, rsize);

   nbad[1] = tst_pdspgvx( n, vecA, mapA, vecB, mapB, vecZ, mapZ, eval,
                          iscratch, isize, iptr, ptr_size ,scratch, rsize);

   nbad[2] = tst_pdspev( n, vecA, mapA, vecZ, mapZ, eval,
                         iscratch, isize, iptr, ptr_size ,scratch, rsize);

   nbad[3] = tst_pdspevx( n, vecA, mapA, vecZ, mapZ, eval,
                          iscratch, isize, iptr, ptr_size ,scratch, rsize);

   mxsync_();
   fflush(stderr);
   fflush(stdout);
   if( me == 0 )
     fprintf( stderr, " \n \n \n" );
   mxsync_();
   fprintf( stderr, " me = %d n = %d pdspgv  nbad = %d \n", me, n, nbad[0]);

   mxsync_();
   if( me == 0 )
     fprintf( stderr, " \n \n \n" );
   mxsync_();
   fprintf( stderr, " me = %d n = %d pdspgvx nbad = %d \n", me, n, nbad[1]);

   mxsync_();
   if( me == 0 )
     fprintf( stderr, " \n \n \n" );
     mxsync_();
   fprintf( stderr, " me = %d n = %d pdspev  nbad = %d \n", me, n, nbad[2]);

   mxsync_();
   if( me == 0 )
     fprintf( stderr, " \n \n \n" );
   mxsync_();
   fprintf( stderr, " me = %d n = %d pdspevx nbad = %d \n", me, n, nbad[3]);

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
  
}

Integer tst_pdspgv( n, vecA, mapA, vecB, mapB, vecZ, mapZ, eval, iscratch,
                    isize, iptr, ptr_size ,scratch, rsize)

  Integer          n, *mapA, *mapB, *mapZ, *iscratch,
                   isize, ptr_size, rsize;
  DoublePrecision  **vecA, **vecB, **vecZ, **iptr;
  DoublePrecision  *eval, *scratch;

{
  Integer          k, ii, indx, nbad, me, ndiff, nprocs, i, info;
  Integer          *inull;
  
  DoublePrecision  *dnull, **ddnull;

  extern Integer mxsync_(), mxmynd_(), mxnprc_();
  extern void    pdspgv();

  me     = mxmynd_();
  nprocs = mxnprc_();

  nbad = 0;

/*--------------------------------------------------
 * Pass pointers to NULL.
 *--------------------------------------------------*/

  
  inull  = NULL;
  dnull  = NULL;
  ddnull = NULL;
  
  indx = 1;  

  fflush(stderr);
  mxsync_();
  if( me == 0 )
    fprintf( stderr, " \n pdspgv 1 = NULL \n");
  pdspgv( inull, &n, vecA, mapA, vecB, mapB, vecZ, mapZ, eval, iscratch,
          &isize, iptr, &ptr_size ,scratch, &rsize, &info);
  if( info != -1 )
    nbad++;
  
  fflush(stderr);
  mxsync_();
  if( me == 0 )
    fprintf( stderr, " \n pdspgv 2 = NULL \n");
  pdspgv( &indx, inull, vecA, mapA, vecB, mapB, vecZ, mapZ, eval, iscratch,
          &isize, iptr, &ptr_size ,scratch, &rsize, &info);
  if( info != -2 )
    nbad++;

  fflush(stderr);
  mxsync_();
  if( me == 0 )
    fprintf( stderr, " \n pdspgv 3 = NULL \n");
  pdspgv( &indx, &n, ddnull, mapA, vecB, mapB, vecZ, mapZ, eval, iscratch,
          &isize, iptr, &ptr_size ,scratch, &rsize, &info);
  if( info != -3 )
    nbad++;

  fflush(stderr);
  mxsync_();
  if( me == 0 )
    fprintf( stderr, " \n pdspgv 4 = NULL \n");
  pdspgv( &indx, &n, vecA, inull, vecB, mapB, vecZ, mapZ, eval, iscratch,
          &isize, iptr, &ptr_size ,scratch, &rsize, &info);
  if( info != -4 )
    nbad++;

  fflush(stderr);
  mxsync_();
  if( me == 0 )
    fprintf( stderr, " \n pdspgv 5 = NULL \n");
  pdspgv( &indx, &n, vecA, mapA, ddnull, mapB, vecZ, mapZ, eval, iscratch,
          &isize, iptr, &ptr_size ,scratch, &rsize, &info);
  if( info != -5 )
    nbad++;

  fflush(stderr);
  mxsync_();
  if( me == 0 )
    fprintf( stderr, " \n pdspgv 6 = NULL \n");
  pdspgv( &indx, &n, vecA, mapA, vecB, inull, vecZ, mapZ, eval, iscratch,
          &isize, iptr, &ptr_size ,scratch, &rsize, &info);
  if( info != -6 )
    nbad++;

  fflush(stderr);
  mxsync_();
  if( me == 0 )
    fprintf( stderr, " \n pdspgv 7 = NULL \n");
  pdspgv( &indx, &n, vecA, mapA, vecB, mapB, ddnull, mapZ, eval, iscratch,
          &isize, iptr, &ptr_size ,scratch, &rsize, &info);
  if( info != -7 )
    nbad++;

  fflush(stderr);
  mxsync_();
  if( me == 0 )
    fprintf( stderr, " \n pdspgv 8 = NULL \n");
  pdspgv( &indx, &n, vecA, mapA, vecB, mapB, vecZ, inull, eval, iscratch,
          &isize, iptr, &ptr_size ,scratch, &rsize, &info);
  if( info != -8 )
    nbad++;

  fflush(stderr);
  mxsync_();
  if( me == 0 )
    fprintf( stderr, " \n pdspgv 9 = NULL \n");
  pdspgv( &indx, &n, vecA, mapA, vecB, mapB, vecZ, mapZ, dnull, iscratch,
          &isize, iptr, &ptr_size ,scratch, &rsize, &info);
  if( info != -9 )
    nbad++;

  fflush(stderr);
  mxsync_();
  if( me == 0 )
    fprintf( stderr, " \n pdspgv 10 = NULL \n");
  pdspgv( &indx, &n, vecA, mapA, vecB, mapB, vecZ, mapZ, eval, inull,
          &isize, iptr, &ptr_size ,scratch, &rsize, &info);
  if( info != -10 )
    nbad++;

  fflush(stderr);
  mxsync_();
  if( me == 0 )
    fprintf( stderr, " \n pdspgv 11 = NULL \n");
  pdspgv( &indx, &n, vecA, mapA, vecB, mapB, vecZ, mapZ, eval, iscratch,
          inull, iptr, &ptr_size ,scratch, &rsize, &info);
  if( info != -11 )
    nbad++;

  fflush(stderr);
  mxsync_();
  if( me == 0 )
    fprintf( stderr, " \n pdspgv 12 = NULL \n");
  pdspgv( &indx, &n, vecA, mapA, vecB, mapB, vecZ, mapZ, eval, iscratch,
          &isize, ddnull, &ptr_size ,scratch, &rsize, &info);
  if( info != -12 )
    nbad++;

  fflush(stderr);
  mxsync_();
  if( me == 0 )
    fprintf( stderr, " \n pdspgv 13 = NULL \n");
  pdspgv( &indx, &n, vecA, mapA, vecB, mapB, vecZ, mapZ, eval, iscratch,
          &isize, iptr, inull,scratch, &rsize, &info);
  if( info != -13 )
    nbad++;

  mxsync_();
  if( me == 0 )
    fprintf( stderr, " \n pdspgv 14 = NULL \n");
  pdspgv( &indx, &n, vecA, mapA, vecB, mapB, vecZ, mapZ, eval, iscratch,
          &isize, iptr, &ptr_size ,dnull, &rsize, &info);
  if( info != -14 )
    nbad++;

  fflush(stderr);
  mxsync_();
  if( me == 0 )
    fprintf( stderr, " \n pdspgv 15 = NULL \n");
  pdspgv( &indx, &n, vecA, mapA, vecB, mapB, vecZ, mapZ, eval, iscratch,
          &isize, iptr, &ptr_size ,scratch, inull, &info);
  if( info != -15 )
    nbad++;

  fflush(stderr);
  mxsync_();
  if( me == 0 )
    fprintf( stderr, " \n pdspgv 16 = NULL \n");
  pdspgv( &indx, &n, vecA, mapA, vecB, mapB, vecZ, mapZ, eval, iscratch,
          &isize, iptr, &ptr_size ,scratch, &rsize, inull );

  /* Can't check info since passed NULL for info. */

  fflush(stderr);
  mxsync_();
  if( me == 0 )
    fprintf( stderr, " \n pdspgv 3 = NULL \n");
  dnull = vecA[0];
  vecA[0] = NULL;
  pdspgv( &indx, &n, vecA, mapA, vecB, mapB, vecZ, mapZ, eval, iscratch,
          &isize, iptr, &ptr_size ,scratch, &rsize, &info );
  if( info != -3 )
    nbad++;
  vecA[0] = dnull;

  mxsync_();
  if( me == 0 )
    fprintf( stderr, " \n pdspgv 5 = NULL \n");
  dnull = vecB[0];
  vecB[0] = NULL;
  pdspgv( &indx, &n, vecA, mapA, vecB, mapB, vecZ, mapZ, eval, iscratch,
          &isize, iptr, &ptr_size ,scratch, &rsize, &info );
  if( info != -5 )
    nbad++;
  vecB[0] = dnull;

  fflush(stderr);
  mxsync_();
  if( me == 0 )
    fprintf( stderr, " \n pdspgv 7 = NULL \n");
  dnull = vecZ[0];
  vecZ[0] = NULL;
  pdspgv( &indx, &n, vecA, mapA, vecB, mapB, vecZ, mapZ, eval, iscratch,
          &isize, iptr, &ptr_size ,scratch, &rsize, &info );
  if( info != -7 )
    nbad++;
  vecZ[0] = dnull;

/*--------------------------------------------------
 * Invalid input
 *--------------------------------------------------*/

  fflush(stderr);
  mxsync_();
  if( me == 0 )
    fprintf( stderr, " \n pdspgv 1 invalid input \n");
  k = 2;
  pdspgv( &k, &n, vecA, mapA, vecB, mapB, vecZ, mapZ, eval, iscratch,
          &isize, iptr, &ptr_size ,scratch, &rsize, &info);
  if( info != -1 )
    nbad++;

  fflush(stderr);
  mxsync_();
  if( me == 0 )
    fprintf( stderr, " \n pdspgv 2 invalid input \n");
  k = -1;
  pdspgv( &indx, &k, vecA, mapA, vecB, mapB, vecZ, mapZ, eval, iscratch,
          &isize, iptr, &ptr_size ,scratch, &rsize, &info);
  if( info != -2 )
    nbad++;

  fflush(stderr);
  mxsync_();
  if( me == 0 )
    fprintf( stderr, " \n pdspgv 4 invalid input \n");
  mapA[0] = -1;
  pdspgv( &indx, &n, vecA, mapA, vecB, mapB, vecZ, mapZ, eval, iscratch,
          &isize, iptr, &ptr_size ,scratch, &rsize, &info);
  if( info != -4 )
    nbad++;
  mapA[0] = 0;

  fflush(stderr);
  mxsync_();
  if( me == 0 )
    fprintf( stderr, " \n pdspgv 6 invalid input \n");
  mapB[0] = -1;
  pdspgv( &indx, &n, vecA, mapA, vecB, mapB, vecZ, mapZ, eval, iscratch,
          &isize, iptr, &ptr_size ,scratch, &rsize, &info);
  if( info != -6 )
    nbad++;
  mapB[0] = 0;

  fflush(stderr);
  mxsync_();
  if( me == 0 )
    fprintf( stderr, " \n pdspgv 8 invalid input \n");
  mapZ[0] = -1;
  pdspgv( &indx, &n, vecA, mapA, vecB, mapB, vecZ, mapZ, eval, iscratch,
          &isize, iptr, &ptr_size ,scratch, &rsize, &info);
  if( info != -8 )
    nbad++;
  mapZ[0] = 0;

  fflush(stderr);
  mxsync_();
  if( me == 0 )
    fprintf( stderr, " \n pdspgv 11 invalid input \n");
  k = 1;
  pdspgv( &indx, &n, vecA, mapA, vecB, mapB, vecZ, mapZ, eval, iscratch,
          &k, iptr, &ptr_size ,scratch, &rsize, &info);
  if( info != -11 )
    nbad++;

  fflush(stderr);
  mxsync_();
  if( me == 0 )
    fprintf( stderr, " \n pdspgv 13 invalid input \n");
  fflush(stderr);
  k = 1;
  pdspgv( &indx, &n, vecA, mapA, vecB, mapB, vecZ, mapZ, eval, iscratch,
          &isize, iptr, &k,scratch, &rsize, &info);
  if( info != -13 )
    nbad++;

  fflush(stderr);
  mxsync_();
  if( me == 0 )
    fprintf( stderr, " \n pdspgv 15 invalid input \n");
  fflush(stderr);
  k = 1;
  pdspgv( &indx, &n, vecA, mapA, vecB, mapB, vecZ, mapZ, eval, iscratch,
          &isize, iptr, &ptr_size ,scratch, &k, &info);
  if( info != -15 )
    nbad++;

  if( nprocs > 1 ) {
    fflush(stderr);
    mxsync_();
    if( me == 0 )
      fprintf( stderr, " \n pdspgv 8 invalid input \n");
    fflush(stderr);
    
    for ( ii = 0;  ii < n; ii++ )
      mapZ[ii] = 0;
    pdspgv( &indx, &n, vecA, mapA, vecB, mapB, vecZ, mapZ, eval, iscratch,
            &isize, iptr, &ptr_size ,scratch, &rsize, &info);
  if( info != -8 )
    nbad++;
    for ( ii = 0;  ii < n; ii++ ) {
      k = ( ii % nprocs);
      mapZ[ii] = k;
    }
  }

/*--------------------------------------------------
 * Data differs on processors
 *--------------------------------------------------*/

  ndiff = nprocs / 2;

  if( me == ndiff && nprocs > 1 ) {

    fflush(stderr);
    mxsync_();
    fprintf( stderr, " \n pdspgv 1 differs on processor = %d  \n", ndiff);
    k = 0;
    pdspgv( &k, &n, vecA, mapA, vecB, mapB, vecZ, mapZ, eval, iscratch,
            &isize, iptr, &ptr_size ,scratch, &rsize, &info);
  if( info != -51 )
    nbad++;

    fflush(stderr);
    mxsync_();
    fprintf( stderr, " \n pdspgv 4 differs on processor = %d  \n", ndiff);
    mapA[0] = 1;
    mapB[0] = 1;
    mapZ[0] = 1;
    mapA[1] = 0;
    mapB[1] = 0;
    mapZ[1] = 0;
    pdspgv( &indx, &n, vecA, mapA, vecB, mapB, vecZ, mapZ, eval, iscratch,
            &isize, iptr, &ptr_size ,scratch, &rsize, &info);
  if( info != -51 )
    nbad++;
    mapA[0] = 0;
    mapB[0] = 0;
    mapZ[0] = 0;
    mapA[1] = 1;
    mapB[1] = 1;
    mapZ[1] = 1;

  } else if( nprocs > 1 ) {

    for ( i = 0; i < 2; i++ ) {
      fflush(stderr);
      mxsync_();
      pdspgv( &indx, &n, vecA, mapA, vecB, mapB, vecZ, mapZ, eval, iscratch,
              &isize, iptr, &ptr_size ,scratch, &rsize, &info);
  if( info != -51 )
    nbad++;
    }
  }
  return( nbad );
}
  
Integer tst_pdspgvx( n, vecA, mapA, vecB, mapB, vecZ, mapZ, eval, iscratch,
                    isize, iptr, ptr_size ,scratch, rsize)

  Integer          n, *mapA, *mapB, *mapZ, *iscratch,
                   isize, ptr_size, rsize;
  DoublePrecision  **vecA, **vecB, **vecZ, **iptr;
  DoublePrecision  *eval, *scratch;

{
  Integer          k, ii, nbad, me, ndiff, nprocs, i, info,
                   ifact, ivector, irange, ilb, iub, meigval;
  Integer          *inull;
  DoublePrecision  lb, ub, abstol, dk;
  
  DoublePrecision  *dnull, **ddnull;

  extern Integer mxsync_(), mxmynd_(), mxnprc_();
  extern void    pdspgvx();

  me     = mxmynd_();
  nprocs = mxnprc_();

  ifact   = 1;
  ivector = 1;
  irange  = 1;
  ilb     = -1;
  iub     = 0;
  meigval = 0;

  lb     = 1.e0;
  ub     = 0.e0;
  abstol = 0.e0;

  nbad = 0;

/*--------------------------------------------------
 * Pass pointers to NULL.
 *--------------------------------------------------*/

  
  inull  = NULL;
  dnull  = NULL;
  ddnull = NULL;
  
  ifact = 1;  

  fflush(stderr);
  mxsync_();
  if( me == 0 )
    fprintf( stderr, " \n pdspgvx 1 = NULL \n");
  pdspgvx( inull, &ivector, &irange, &n, vecA, mapA, vecB, mapB,
           &lb, &ub, &ilb, &iub, &abstol, &meigval, vecZ, mapZ, eval, iscratch,
          &isize, iptr, &ptr_size ,scratch, &rsize, &info);
  if( info != -1 )
    nbad++;
  
  fflush(stderr);
  mxsync_();
  if( me == 0 )
    fprintf( stderr, " \n pdspgvx 2 = NULL \n");
  pdspgvx( &ifact, inull, &irange, &n, vecA, mapA, vecB, mapB,
           &lb, &ub, &ilb, &iub, &abstol, &meigval, vecZ, mapZ, eval, iscratch,
          &isize, iptr, &ptr_size ,scratch, &rsize, &info);
  if( info != -2 )
    nbad++;

  fflush(stderr);
  mxsync_();
  if( me == 0 )
    fprintf( stderr, " \n pdspgvx 3 = NULL \n");
  pdspgvx( &ifact, &ivector, inull, &n, vecA, mapA, vecB, mapB,
           &lb, &ub, &ilb, &iub, &abstol, &meigval, vecZ, mapZ, eval, iscratch,
          &isize, iptr, &ptr_size ,scratch, &rsize, &info);
  if( info != -3 )
    nbad++;

  fflush(stderr);
  mxsync_();
  if( me == 0 )
    fprintf( stderr, " \n pdspgvx 4 = NULL \n");
  pdspgvx( &ifact, &ivector, &irange, inull, vecA, mapA, vecB, mapB,
           &lb, &ub, &ilb, &iub, &abstol, &meigval, vecZ, mapZ, eval, iscratch,
          &isize, iptr, &ptr_size ,scratch, &rsize, &info);
  if( info != -4 )
    nbad++;

  fflush(stderr);
  mxsync_();
  if( me == 0 )
    fprintf( stderr, " \n pdspgvx 5 = NULL \n");
  pdspgvx( &ifact, &ivector, &irange, &n, ddnull, mapA, vecB, mapB,
           &lb, &ub, &ilb, &iub, &abstol, &meigval, vecZ, mapZ, eval, iscratch,
          &isize, iptr, &ptr_size ,scratch, &rsize, &info);
  if( info != -5 )
    nbad++;

  fflush(stderr);
  mxsync_();
  if( me == 0 )
    fprintf( stderr, " \n pdspgvx 6 = NULL \n");
  pdspgvx( &ifact, &ivector, &irange, &n, vecA, inull, vecB, mapB,
           &lb, &ub, &ilb, &iub, &abstol, &meigval, vecZ, mapZ, eval, iscratch,
          &isize, iptr, &ptr_size ,scratch, &rsize, &info);
  if( info != -6 )
    nbad++;

  fflush(stderr);
  mxsync_();
  if( me == 0 )
    fprintf( stderr, " \n pdspgvx 7 = NULL \n");
  pdspgvx( &ifact, &ivector, &irange, &n, vecA, mapA, ddnull, mapB,
           &lb, &ub, &ilb, &iub, &abstol, &meigval, vecZ, mapZ, eval, iscratch,
          &isize, iptr, &ptr_size ,scratch, &rsize, &info);
  if( info != -7 )
    nbad++;

  fflush(stderr);
  mxsync_();
  if( me == 0 )
    fprintf( stderr, " \n pdspgvx 8 = NULL \n");
  pdspgvx( &ifact, &ivector, &irange, &n, vecA, mapA, vecB, inull,
           &lb, &ub, &ilb, &iub, &abstol, &meigval, vecZ, mapZ, eval, iscratch,
          &isize, iptr, &ptr_size ,scratch, &rsize, &info);
  if( info != -8 )
    nbad++;

  fflush(stderr);
  mxsync_();
  if( me == 0 )
    fprintf( stderr, " \n pdspgvx 9 = NULL \n");
  pdspgvx( &ifact, &ivector, &irange, &n, vecA, mapA, vecB, mapB,
           dnull, &ub, &ilb, &iub, &abstol, &meigval, vecZ, mapZ, eval, iscratch,
          &isize, iptr, &ptr_size ,scratch, &rsize, &info);
  if( info != -9 )
    nbad++;

  fflush(stderr);
  mxsync_();
  if( me == 0 )
    fprintf( stderr, " \n pdspgvx 10 = NULL \n");
  pdspgvx( &ifact, &ivector, &irange, &n, vecA, mapA, vecB, mapB,
           &lb, dnull, &ilb, &iub, &abstol, &meigval, vecZ, mapZ, eval, iscratch,
          &isize, iptr, &ptr_size ,scratch, &rsize, &info);
  if( info != -10 )
    nbad++;

  fflush(stderr);
  mxsync_();
  if( me == 0 )
    fprintf( stderr, " \n pdspgvx 11 = NULL \n");
  pdspgvx( &ifact, &ivector, &irange, &n, vecA, mapA, vecB, mapB,
           &lb, &ub, inull, &iub, &abstol, &meigval, vecZ, mapZ, eval, iscratch,
          &isize, iptr, &ptr_size ,scratch, &rsize, &info);
  if( info != -11 )
    nbad++;

  fflush(stderr);
  mxsync_();
  if( me == 0 )
    fprintf( stderr, " \n pdspgvx 12 = NULL \n");
  pdspgvx( &ifact, &ivector, &irange, &n, vecA, mapA, vecB, mapB,
           &lb, &ub, &ilb, inull, &abstol, &meigval, vecZ, mapZ, eval, iscratch,
          &isize, iptr, &ptr_size ,scratch, &rsize, &info);
  if( info != -12 )
    nbad++;

  fflush(stderr);
  mxsync_();
  if( me == 0 )
    fprintf( stderr, " \n pdspgvx 13 = NULL \n");
  pdspgvx( &ifact, &ivector, &irange, &n, vecA, mapA, vecB, mapB,
           &lb, &ub, &ilb, &iub, dnull, &meigval, vecZ, mapZ, eval, iscratch,
          &isize, iptr, &ptr_size ,scratch, &rsize, &info);
  if( info != -13 )
    nbad++;

  fflush(stderr);
  mxsync_();
  if( me == 0 )
    fprintf( stderr, " \n pdspgvx 14 = NULL \n");
  pdspgvx( &ifact, &ivector, &irange, &n, vecA, mapA, vecB, mapB,
           &lb, &ub, &ilb, &iub, &abstol, inull, vecZ, mapZ, eval, iscratch,
          &isize, iptr, &ptr_size ,scratch, &rsize, &info);
  if( info != -14 )
    nbad++;

  fflush(stderr);
  mxsync_();
  if( me == 0 )
    fprintf( stderr, " \n pdspgvx 15 = NULL \n");
  pdspgvx( &ifact, &ivector, &irange, &n, vecA, mapA, vecB, mapB,
           &lb, &ub, &ilb, &iub, &abstol, &meigval, ddnull, mapZ, eval, iscratch,
          &isize, iptr, &ptr_size ,scratch, &rsize, &info);
  if( info != -15 )
    nbad++;

  fflush(stderr);
  mxsync_();
  if( me == 0 )
    fprintf( stderr, " \n pdspgvx 16 = NULL \n");
  pdspgvx( &ifact, &ivector, &irange, &n, vecA, mapA, vecB, mapB,
           &lb, &ub, &ilb, &iub, &abstol, &meigval, vecZ, inull, eval, iscratch,
          &isize, iptr, &ptr_size ,scratch, &rsize, &info);
  if( info != -16 )
    nbad++;

  fflush(stderr);
  mxsync_();
  if( me == 0 )
    fprintf( stderr, " \n pdspgvx 17 = NULL \n");
  pdspgvx( &ifact, &ivector, &irange, &n, vecA, mapA, vecB, mapB,
           &lb, &ub, &ilb, &iub, &abstol, &meigval, vecZ, mapZ, dnull, iscratch,
          &isize, iptr, &ptr_size ,scratch, &rsize, &info);
  if( info != -17 )
    nbad++;

  fflush(stderr);
  mxsync_();
  if( me == 0 )
    fprintf( stderr, " \n pdspgvx 18 = NULL \n");
  pdspgvx( &ifact, &ivector, &irange, &n, vecA, mapA, vecB, mapB,
           &lb, &ub, &ilb, &iub, &abstol, &meigval, vecZ, mapZ, eval, inull,
          &isize, iptr, &ptr_size ,scratch, &rsize, &info);
  if( info != -18 )
    nbad++;

  fflush(stderr);
  mxsync_();
  if( me == 0 )
    fprintf( stderr, " \n pdspgvx 19 = NULL \n");
  pdspgvx( &ifact, &ivector, &irange, &n, vecA, mapA, vecB, mapB,
           &lb, &ub, &ilb, &iub, &abstol, &meigval, vecZ, mapZ, eval, iscratch,
          inull, iptr, &ptr_size ,scratch, &rsize, &info);
  if( info != -19 )
    nbad++;

  fflush(stderr);
  mxsync_();
  if( me == 0 )
    fprintf( stderr, " \n pdspgvx 20 = NULL \n");
  pdspgvx( &ifact, &ivector, &irange, &n, vecA, mapA, vecB, mapB,
           &lb, &ub, &ilb, &iub, &abstol, &meigval, vecZ, mapZ, eval, iscratch,
          &isize, ddnull, &ptr_size ,scratch, &rsize, &info);
  if( info != -20 )
    nbad++;

  fflush(stderr);
  mxsync_();
  if( me == 0 )
    fprintf( stderr, " \n pdspgvx 21 = NULL \n");
  pdspgvx( &ifact, &ivector, &irange, &n, vecA, mapA, vecB, mapB,
           &lb, &ub, &ilb, &iub, &abstol, &meigval, vecZ, mapZ, eval, iscratch,
          &isize, iptr, inull,scratch, &rsize, &info);
  if( info != -21 )
    nbad++;

  mxsync_();
  if( me == 0 )
    fprintf( stderr, " \n pdspgvx 22 = NULL \n");
  pdspgvx( &ifact, &ivector, &irange, &n, vecA, mapA, vecB, mapB,
           &lb, &ub, &ilb, &iub, &abstol, &meigval, vecZ, mapZ, eval, iscratch,
          &isize, iptr, &ptr_size ,dnull, &rsize, &info);
  if( info != -22 )
    nbad++;

  fflush(stderr);
  mxsync_();
  if( me == 0 )
    fprintf( stderr, " \n pdspgvx 23 = NULL \n");
  pdspgvx( &ifact, &ivector, &irange, &n, vecA, mapA, vecB, mapB,
           &lb, &ub, &ilb, &iub, &abstol, &meigval, vecZ, mapZ, eval, iscratch,
          &isize, iptr, &ptr_size ,scratch, inull, &info);
  if( info != -23 )
    nbad++;

  fflush(stderr);
  mxsync_();
  if( me == 0 )
    fprintf( stderr, " \n pdspgvx 24 = NULL \n");
  pdspgvx( &ifact, &ivector, &irange, &n, vecA, mapA, vecB, mapB,
           &lb, &ub, &ilb, &iub, &abstol, &meigval, vecZ, mapZ, eval, iscratch,
          &isize, iptr, &ptr_size ,scratch, &rsize, inull );

  /* Can't check info since passed NULL for info. */

  fflush(stderr);
  mxsync_();
  if( me == 0 )
    fprintf( stderr, " \n pdspgvx 5 = NULL \n");
  dnull = vecA[0];
  vecA[0] = NULL;
  pdspgvx( &ifact, &ivector, &irange, &n, vecA, mapA, vecB, mapB,
           &lb, &ub, &ilb, &iub, &abstol, &meigval, vecZ, mapZ, eval, iscratch,
          &isize, iptr, &ptr_size ,scratch, &rsize, &info );
  if( info != -5 )
    nbad++;
  vecA[0] = dnull;

  mxsync_();
  if( me == 0 )
    fprintf( stderr, " \n pdspgvx 7 = NULL \n");
  dnull = vecB[0];
  vecB[0] = NULL;
  pdspgvx( &ifact, &ivector, &irange, &n, vecA, mapA, vecB, mapB,
           &lb, &ub, &ilb, &iub, &abstol, &meigval, vecZ, mapZ, eval, iscratch,
          &isize, iptr, &ptr_size ,scratch, &rsize, &info );
  if( info != -7 )
    nbad++;
  vecB[0] = dnull;

  fflush(stderr);
  mxsync_();
  if( me == 0 )
    fprintf( stderr, " \n pdspgvx 15 = NULL \n");
  dnull = vecZ[0];
  vecZ[0] = NULL;
  pdspgvx( &ifact, &ivector, &irange, &n, vecA, mapA, vecB, mapB,
           &lb, &ub, &ilb, &iub, &abstol, &meigval, vecZ, mapZ, eval, iscratch,
          &isize, iptr, &ptr_size ,scratch, &rsize, &info );
  if( info != -15 )
    nbad++;
  vecZ[0] = dnull;

/*--------------------------------------------------
 * Invalid input
 *--------------------------------------------------*/

  fflush(stderr);
  mxsync_();
  if( me == 0 )
    fprintf( stderr, " \n pdspgvx 1 invalid input \n");
  k = 2;
  pdspgvx( &k, &ivector, &irange, &n, vecA, mapA, vecB, mapB,
           &lb, &ub, &ilb, &iub, &abstol, &meigval, vecZ, mapZ, eval, iscratch,
          &isize, iptr, &ptr_size ,scratch, &rsize, &info);
  if( info != -1 )
    nbad++;

  fflush(stderr);
  mxsync_();
  if( me == 0 )
    fprintf( stderr, " \n pdspgvx 2 invalid input \n");
  k = -1;
  pdspgvx( &ifact, &k, &irange, &n, vecA, mapA, vecB, mapB,
           &lb, &ub, &ilb, &iub, &abstol, &meigval, vecZ, mapZ, eval, iscratch,
          &isize, iptr, &ptr_size ,scratch, &rsize, &info);
  if( info != -2 )
    nbad++;

  fflush(stderr);
  mxsync_();
  if( me == 0 )
    fprintf( stderr, " \n pdspgvx 3 invalid input \n");
  k = -1;
  pdspgvx( &ifact, &ivector, &k, &n, vecA, mapA, vecB, mapB,
           &lb, &ub, &ilb, &iub, &abstol, &meigval, vecZ, mapZ, eval, iscratch,
          &isize, iptr, &ptr_size ,scratch, &rsize, &info);
  if( info != -3 )
    nbad++;

  fflush(stderr);
  mxsync_();
  if( me == 0 )
    fprintf( stderr, " \n pdspgvx 4 invalid input \n");
  k = -1;
  pdspgvx( &ifact, &ivector, &irange, &k, vecA, mapA, vecB, mapB,
           &lb, &ub, &ilb, &iub, &abstol, &meigval, vecZ, mapZ, eval, iscratch,
          &isize, iptr, &ptr_size ,scratch, &rsize, &info);
  if( info != -4 )
    nbad++;

  fflush(stderr);
  mxsync_();
  if( me == 0 )
    fprintf( stderr, " \n pdspgvx 6 invalid input \n");
  mapA[0] = -1;
  pdspgvx( &ifact, &ivector, &irange, &n, vecA, mapA, vecB, mapB,
           &lb, &ub, &ilb, &iub, &abstol, &meigval, vecZ, mapZ, eval, iscratch,
          &isize, iptr, &ptr_size ,scratch, &rsize, &info);
  if( info != -6 )
    nbad++;
  mapA[0] = 0;

  fflush(stderr);
  mxsync_();
  if( me == 0 )
    fprintf( stderr, " \n pdspgvx 8 invalid input \n");
  mapB[0] = -1;
  pdspgvx( &ifact, &ivector, &irange, &n, vecA, mapA, vecB, mapB,
           &lb, &ub, &ilb, &iub, &abstol, &meigval, vecZ, mapZ, eval, iscratch,
          &isize, iptr, &ptr_size ,scratch, &rsize, &info);
  if( info != -8 )
    nbad++;
  mapB[0] = 0;

  fflush(stderr);
  mxsync_();
  if( me == 0 )
    fprintf( stderr, " \n pdspgvx 10 invalid input \n");
  k = 2;
  pdspgvx( &ifact, &ivector, &k, &n, vecA, mapA, vecB, mapB,
           &lb, &ub, &ilb, &iub, &abstol, &meigval, vecZ, mapZ, eval, iscratch,
          &isize, iptr, &ptr_size ,scratch, &rsize, &info);
  if( info != -10 )
    nbad++;

  fflush(stderr);
  mxsync_();
  if( me == 0 )
    fprintf( stderr, " \n pdspgvx 11 invalid input \n");
  k = 3;
  ilb = -1;
  iub = 0;
  pdspgvx( &ifact, &ivector, &k, &n, vecA, mapA, vecB, mapB,
           &lb, &ub, &ilb, &iub, &abstol, &meigval, vecZ, mapZ, eval, iscratch,
          &isize, iptr, &ptr_size ,scratch, &rsize, &info);
  if( info != -11 )
    nbad++;

  fflush(stderr);
  mxsync_();
  if( me == 0 )
    fprintf( stderr, " \n pdspgvx 12 invalid input \n");
  k = 3;
  ilb = 2;
  iub = 1;
  pdspgvx( &ifact, &ivector, &k, &n, vecA, mapA, vecB, mapB,
           &lb, &ub, &ilb, &iub, &abstol, &meigval, vecZ, mapZ, eval, iscratch,
          &isize, iptr, &ptr_size ,scratch, &rsize, &info);
  if( info != -12 )
    nbad++;

  fflush(stderr);
  mxsync_();
  if( me == 0 )
    fprintf( stderr, " \n pdspgvx 12 invalid input \n");
  k = 3;
  ilb = 1;
  iub = n+1;
  pdspgvx( &ifact, &ivector, &k, &n, vecA, mapA, vecB, mapB,
           &lb, &ub, &ilb, &iub, &abstol, &meigval, vecZ, mapZ, eval, iscratch,
          &isize, iptr, &ptr_size ,scratch, &rsize, &info);
  if( info != -12 )
    nbad++;

  fflush(stderr);
  mxsync_();
  if( me == 0 )
    fprintf( stderr, " \n pdspgvx 16 invalid input \n");
  mapZ[0] = -1;
  pdspgvx( &ifact, &ivector, &irange, &n, vecA, mapA, vecB, mapB,
           &lb, &ub, &ilb, &iub, &abstol, &meigval, vecZ, mapZ, eval, iscratch,
          &isize, iptr, &ptr_size ,scratch, &rsize, &info);
  if( info != -16 )
    nbad++;
  mapZ[0] = 0;

  fflush(stderr);
  mxsync_();
  if( me == 0 )
    fprintf( stderr, " \n pdspgvx 19 invalid input \n");
  k = 1;
  pdspgvx( &ifact, &ivector, &irange, &n, vecA, mapA, vecB, mapB,
           &lb, &ub, &ilb, &iub, &abstol, &meigval, vecZ, mapZ, eval, iscratch,
          &k, iptr, &ptr_size ,scratch, &rsize, &info);
  if( info != -19 )
    nbad++;

  fflush(stderr);
  mxsync_();
  if( me == 0 )
    fprintf( stderr, " \n pdspgvx 21 invalid input \n");
  fflush(stderr);
  k = 1;
  pdspgvx( &ifact, &ivector, &irange, &n, vecA, mapA, vecB, mapB,
           &lb, &ub, &ilb, &iub, &abstol, &meigval, vecZ, mapZ, eval, iscratch,
          &isize, iptr, &k,scratch, &rsize, &info);
  if( info != -21 )
    nbad++;

  fflush(stderr);
  mxsync_();
  if( me == 0 )
    fprintf( stderr, " \n pdspgvx 23 invalid input \n");
  fflush(stderr);
  k = 1;
  pdspgvx( &ifact, &ivector, &irange, &n, vecA, mapA, vecB, mapB,
           &lb, &ub, &ilb, &iub, &abstol, &meigval, vecZ, mapZ, eval, iscratch,
          &isize, iptr, &ptr_size ,scratch, &k, &info);
  if( info != -23 )
    nbad++;

  if( nprocs > 1 ) {
    fflush(stderr);
    mxsync_();
    if( me == 0 )
      fprintf( stderr, " \n pdspgvx 16 invalid input \n");
    fflush(stderr);
    
    for ( ii = 0;  ii < n; ii++ )
      mapZ[ii] = 0;
    pdspgvx( &ifact, &ivector, &irange, &n, vecA, mapA, vecB, mapB,
           &lb, &ub, &ilb, &iub, &abstol, &meigval, vecZ, mapZ, eval, iscratch,
            &isize, iptr, &ptr_size ,scratch, &rsize, &info);
  if( info != -16 )
    nbad++;
    for ( ii = 0;  ii < n; ii++ ) {
      k = ( ii % nprocs);
      mapZ[ii] = k;
    }
  }

/*--------------------------------------------------
 * Data differs on processors
 *--------------------------------------------------*/

  ndiff = nprocs / 2;

  if( me == ndiff && nprocs > 1 ) {

    fflush(stderr);
    mxsync_();
    fprintf( stderr, " \n pdspgvx 1 differs on processor = %d  \n", ndiff);
    k = 0;
    pdspgvx( &k, &ivector, &irange, &n, vecA, mapA, vecB, mapB, 
          &lb, &ub, &ilb, &iub, &abstol, &meigval, vecZ, mapZ, eval, iscratch,
            &isize, iptr, &ptr_size ,scratch, &rsize, &info);
    if( info != -51 )
      nbad++;

    fflush(stderr);
    mxsync_();
    fprintf( stderr, " \n pdspgvx 9 on processor = %d  \n", ndiff);
    dk = -10.e0;
    pdspgvx( &ifact, &ivector, &irange, &n, vecA, mapA, vecB, mapB, 
          &dk, &ub, &ilb, &iub, &abstol, &meigval, vecZ, mapZ, eval, iscratch,
            &isize, iptr, &ptr_size ,scratch, &rsize, &info);
    if( info != -51 )
      nbad++;

    fflush(stderr);
    mxsync_();
    fprintf( stderr, " \n pdspgvx maps differ on processor = %d  \n", ndiff);
    mapA[0] = 1;
    mapB[0] = 1;
    mapZ[0] = 1;
    mapA[1] = 0;
    mapB[1] = 0;
    mapZ[1] = 0;
    pdspgvx( &ifact, &ivector, &irange, &n, vecA, mapA, vecB, mapB,
           &lb, &ub, &ilb, &iub, &abstol, &meigval, vecZ, mapZ, eval, iscratch,
            &isize, iptr, &ptr_size ,scratch, &rsize, &info);
  if( info != -51 )
    nbad++;
    mapA[0] = 0;
    mapB[0] = 0;
    mapZ[0] = 0;
    mapA[1] = 1;
    mapB[1] = 1;
    mapZ[1] = 1;

  } else if( nprocs > 1 ) {

    for ( i = 0; i < 3; i++ ) {
      fflush(stderr);
      mxsync_();
      pdspgvx( &ifact, &ivector, &irange, &n, vecA, mapA, vecB, mapB,
           &lb, &ub, &ilb, &iub, &abstol, &meigval, vecZ, mapZ, eval, iscratch,
              &isize, iptr, &ptr_size ,scratch, &rsize, &info);
  if( info != -51 )
    nbad++;
    }
  }
  return( nbad );
}
Integer tst_pdspevx( n, vecA, mapA,  vecZ, mapZ, eval, iscratch,
                    isize, iptr, ptr_size ,scratch, rsize)

  Integer          n, *mapA, *mapZ, *iscratch,
                   isize, ptr_size, rsize;
  DoublePrecision  **vecA, **vecZ, **iptr;
  DoublePrecision  *eval, *scratch;

{
  Integer          k, ii, nbad, me, ndiff, nprocs, i, info,
                   ivector, irange, ilb, iub, meigval;
  Integer          *inull;
  DoublePrecision  lb, ub, abstol, dk;
  
  DoublePrecision  *dnull, **ddnull;

  extern Integer mxsync_(), mxmynd_(), mxnprc_();
  extern void    pdspevx();

  me     = mxmynd_();
  nprocs = mxnprc_();

  ivector = 1;
  irange  = 1;
  ilb     = -1;
  iub     = 0;
  meigval = 0;

  lb     = 1.e0;
  ub     = 0.e0;
  abstol = 0.e0;

  nbad = 0;

/*--------------------------------------------------
 * Pass pointers to NULL.
 *--------------------------------------------------*/

  
  inull  = NULL;
  dnull  = NULL;
  ddnull = NULL;
  
  fflush(stderr);
  mxsync_();
  if( me == 0 )
    fprintf( stderr, " \n pdspevx 1 = NULL \n");
  pdspevx(  inull, &irange, &n, vecA, mapA, 
           &lb, &ub, &ilb, &iub, &abstol, &meigval, vecZ, mapZ, eval, iscratch,
          &isize, iptr, &ptr_size ,scratch, &rsize, &info);
  if( info != -1 )
    nbad++;

  fflush(stderr);
  mxsync_();
  if( me == 0 )
    fprintf( stderr, " \n pdspevx 2 = NULL \n");
  pdspevx(  &ivector, inull, &n, vecA, mapA, 
           &lb, &ub, &ilb, &iub, &abstol, &meigval, vecZ, mapZ, eval, iscratch,
          &isize, iptr, &ptr_size ,scratch, &rsize, &info);
  if( info != -2 )
    nbad++;

  fflush(stderr);
  mxsync_();
  if( me == 0 )
    fprintf( stderr, " \n pdspevx 3 = NULL \n");
  pdspevx(  &ivector, &irange, inull, vecA, mapA, 
           &lb, &ub, &ilb, &iub, &abstol, &meigval, vecZ, mapZ, eval, iscratch,
          &isize, iptr, &ptr_size ,scratch, &rsize, &info);
  if( info != -3 )
    nbad++;

  fflush(stderr);
  mxsync_();
  if( me == 0 )
    fprintf( stderr, " \n pdspevx 4 = NULL \n");
  pdspevx(  &ivector, &irange, &n, ddnull, mapA, 
           &lb, &ub, &ilb, &iub, &abstol, &meigval, vecZ, mapZ, eval, iscratch,
          &isize, iptr, &ptr_size ,scratch, &rsize, &info);
  if( info != -4 )
    nbad++;

  fflush(stderr);
  mxsync_();
  if( me == 0 )
    fprintf( stderr, " \n pdspevx 5 = NULL \n");
  pdspevx(  &ivector, &irange, &n, vecA, inull, 
           &lb, &ub, &ilb, &iub, &abstol, &meigval, vecZ, mapZ, eval, iscratch,
          &isize, iptr, &ptr_size ,scratch, &rsize, &info);
  if( info != -5 )
    nbad++;

  fflush(stderr);
  mxsync_();
  if( me == 0 )
    fprintf( stderr, " \n pdspevx 6 = NULL \n");
  pdspevx(  &ivector, &irange, &n, vecA, mapA, 
           dnull, &ub, &ilb, &iub, &abstol, &meigval, vecZ, mapZ, eval, iscratch,
          &isize, iptr, &ptr_size ,scratch, &rsize, &info);
  if( info != -6 )
    nbad++;

  fflush(stderr);
  mxsync_();
  if( me == 0 )
    fprintf( stderr, " \n pdspevx 7  = NULL \n");
  pdspevx(  &ivector, &irange, &n, vecA, mapA, 
           &lb, dnull, &ilb, &iub, &abstol, &meigval, vecZ, mapZ, eval, iscratch,
          &isize, iptr, &ptr_size ,scratch, &rsize, &info);
  if( info != -7  )
    nbad++;

  fflush(stderr);
  mxsync_();
  if( me == 0 )
    fprintf( stderr, " \n pdspevx 8  = NULL \n");
  pdspevx(  &ivector, &irange, &n, vecA, mapA, 
           &lb, &ub, inull, &iub, &abstol, &meigval, vecZ, mapZ, eval, iscratch,
          &isize, iptr, &ptr_size ,scratch, &rsize, &info);
  if( info != -8  )
    nbad++;

  fflush(stderr);
  mxsync_();
  if( me == 0 )
    fprintf( stderr, " \n pdspevx 9  = NULL \n");
  pdspevx(  &ivector, &irange, &n, vecA, mapA, 
           &lb, &ub, &ilb, inull, &abstol, &meigval, vecZ, mapZ, eval, iscratch,
          &isize, iptr, &ptr_size ,scratch, &rsize, &info);
  if( info != -9  )
    nbad++;

  fflush(stderr);
  mxsync_();
  if( me == 0 )
    fprintf( stderr, " \n pdspevx 10 = NULL \n");
  pdspevx(  &ivector, &irange, &n, vecA, mapA, 
           &lb, &ub, &ilb, &iub, dnull, &meigval, vecZ, mapZ, eval, iscratch,
          &isize, iptr, &ptr_size ,scratch, &rsize, &info);
  if( info != -10 )
    nbad++;

  fflush(stderr);
  mxsync_();
  if( me == 0 )
    fprintf( stderr, " \n pdspevx 11 = NULL \n");
  pdspevx(  &ivector, &irange, &n, vecA, mapA, 
           &lb, &ub, &ilb, &iub, &abstol, inull, vecZ, mapZ, eval, iscratch,
          &isize, iptr, &ptr_size ,scratch, &rsize, &info);
  if( info != -11 )
    nbad++;

  fflush(stderr);
  mxsync_();
  if( me == 0 )
    fprintf( stderr, " \n pdspevx 12 = NULL \n");
  pdspevx(  &ivector, &irange, &n, vecA, mapA, 
           &lb, &ub, &ilb, &iub, &abstol, &meigval, ddnull, mapZ, eval, iscratch,
          &isize, iptr, &ptr_size ,scratch, &rsize, &info);
  if( info != -12 )
    nbad++;

  fflush(stderr);
  mxsync_();
  if( me == 0 )
    fprintf( stderr, " \n pdspevx 13 = NULL \n");
  pdspevx(  &ivector, &irange, &n, vecA, mapA, 
           &lb, &ub, &ilb, &iub, &abstol, &meigval, vecZ, inull, eval, iscratch,
          &isize, iptr, &ptr_size ,scratch, &rsize, &info);
  if( info != -13 )
    nbad++;

  fflush(stderr);
  mxsync_();
  if( me == 0 )
    fprintf( stderr, " \n pdspevx 14 = NULL \n");
  pdspevx(  &ivector, &irange, &n, vecA, mapA, 
           &lb, &ub, &ilb, &iub, &abstol, &meigval, vecZ, mapZ, dnull, iscratch,
          &isize, iptr, &ptr_size ,scratch, &rsize, &info);
  if( info != -14 )
    nbad++;

  fflush(stderr);
  mxsync_();
  if( me == 0 )
    fprintf( stderr, " \n pdspevx 15 = NULL \n");
  pdspevx(  &ivector, &irange, &n, vecA, mapA, 
           &lb, &ub, &ilb, &iub, &abstol, &meigval, vecZ, mapZ, eval, inull,
          &isize, iptr, &ptr_size ,scratch, &rsize, &info);
  if( info != -15 )
    nbad++;

  fflush(stderr);
  mxsync_();
  if( me == 0 )
    fprintf( stderr, " \n pdspevx 16 = NULL \n");
  pdspevx(  &ivector, &irange, &n, vecA, mapA, 
           &lb, &ub, &ilb, &iub, &abstol, &meigval, vecZ, mapZ, eval, iscratch,
          inull, iptr, &ptr_size ,scratch, &rsize, &info);
  if( info != -16 )
    nbad++;

  fflush(stderr);
  mxsync_();
  if( me == 0 )
    fprintf( stderr, " \n pdspevx 17 = NULL \n");
  pdspevx(  &ivector, &irange, &n, vecA, mapA, 
           &lb, &ub, &ilb, &iub, &abstol, &meigval, vecZ, mapZ, eval, iscratch,
          &isize, ddnull, &ptr_size ,scratch, &rsize, &info);
  if( info != -17 )
    nbad++;

  fflush(stderr);
  mxsync_();
  if( me == 0 )
    fprintf( stderr, " \n pdspevx 18 = NULL \n");
  pdspevx(  &ivector, &irange, &n, vecA, mapA, 
           &lb, &ub, &ilb, &iub, &abstol, &meigval, vecZ, mapZ, eval, iscratch,
          &isize, iptr, inull,scratch, &rsize, &info);
  if( info != -18 )
    nbad++;

  mxsync_();
  if( me == 0 )
    fprintf( stderr, " \n pdspevx 19 = NULL \n");
  pdspevx(  &ivector, &irange, &n, vecA, mapA, 
           &lb, &ub, &ilb, &iub, &abstol, &meigval, vecZ, mapZ, eval, iscratch,
          &isize, iptr, &ptr_size ,dnull, &rsize, &info);
  if( info != -19 )
    nbad++;

  fflush(stderr);
  mxsync_();
  if( me == 0 )
    fprintf( stderr, " \n pdspevx 20 = NULL \n");
  pdspevx(  &ivector, &irange, &n, vecA, mapA, 
           &lb, &ub, &ilb, &iub, &abstol, &meigval, vecZ, mapZ, eval, iscratch,
          &isize, iptr, &ptr_size ,scratch, inull, &info);
  if( info != -20 )
    nbad++;

  fflush(stderr);
  mxsync_();
  if( me == 0 )
    fprintf( stderr, " \n pdspevx 21 = NULL \n");
  pdspevx(  &ivector, &irange, &n, vecA, mapA, 
           &lb, &ub, &ilb, &iub, &abstol, &meigval, vecZ, mapZ, eval, iscratch,
          &isize, iptr, &ptr_size ,scratch, &rsize, inull );

  /* Can't check info since passed NULL for info. */

  fflush(stderr);
  mxsync_();
  if( me == 0 )
    fprintf( stderr, " \n pdspevx 4 = NULL \n");
  dnull = vecA[0];
  vecA[0] = NULL;
  pdspevx(  &ivector, &irange, &n, vecA, mapA, 
           &lb, &ub, &ilb, &iub, &abstol, &meigval, vecZ, mapZ, eval, iscratch,
          &isize, iptr, &ptr_size ,scratch, &rsize, &info );
  if( info != -4 )
    nbad++;
  vecA[0] = dnull;

  fflush(stderr);
  mxsync_();
  if( me == 0 )
    fprintf( stderr, " \n pdspevx 12 = NULL \n");
  dnull = vecZ[0];
  vecZ[0] = NULL;
  pdspevx(  &ivector, &irange, &n, vecA, mapA, 
           &lb, &ub, &ilb, &iub, &abstol, &meigval, vecZ, mapZ, eval, iscratch,
          &isize, iptr, &ptr_size ,scratch, &rsize, &info );
  if( info != -12 )
    nbad++;
  vecZ[0] = dnull;

/*--------------------------------------------------
 * Invalid input
 *--------------------------------------------------*/

  fflush(stderr);
  mxsync_();
  if( me == 0 )
    fprintf( stderr, " \n pdspevx 1 invalid input \n");
  k = -1;
  pdspevx(  &k, &irange, &n, vecA, mapA, 
           &lb, &ub, &ilb, &iub, &abstol, &meigval, vecZ, mapZ, eval, iscratch,
          &isize, iptr, &ptr_size ,scratch, &rsize, &info);
  if( info != -1 )
    nbad++;

  fflush(stderr);
  mxsync_();
  if( me == 0 )
    fprintf( stderr, " \n pdspevx 2 invalid input \n");
  k = -1;
  pdspevx(  &ivector, &k, &n, vecA, mapA, 
           &lb, &ub, &ilb, &iub, &abstol, &meigval, vecZ, mapZ, eval, iscratch,
          &isize, iptr, &ptr_size ,scratch, &rsize, &info);
  if( info != -2 )
    nbad++;

  fflush(stderr);
  mxsync_();
  if( me == 0 )
    fprintf( stderr, " \n pdspevx 3 invalid input \n");
  k = -1;
  pdspevx(  &ivector, &irange, &k, vecA, mapA, 
           &lb, &ub, &ilb, &iub, &abstol, &meigval, vecZ, mapZ, eval, iscratch,
          &isize, iptr, &ptr_size ,scratch, &rsize, &info);
  if( info != -3 )
    nbad++;

  fflush(stderr);
  mxsync_();
  if( me == 0 )
    fprintf( stderr, " \n pdspevx 5 invalid input \n");
  mapA[0] = -1;
  pdspevx(  &ivector, &irange, &n, vecA, mapA, 
           &lb, &ub, &ilb, &iub, &abstol, &meigval, vecZ, mapZ, eval, iscratch,
          &isize, iptr, &ptr_size ,scratch, &rsize, &info);
  if( info != -5 )
    nbad++;
  mapA[0] = 0;

  fflush(stderr);
  mxsync_();
  if( me == 0 )
    fprintf( stderr, " \n pdspevx 7  invalid input \n");
  k = 2;
  pdspevx(  &ivector, &k, &n, vecA, mapA, 
           &lb, &ub, &ilb, &iub, &abstol, &meigval, vecZ, mapZ, eval, iscratch,
          &isize, iptr, &ptr_size ,scratch, &rsize, &info);
  if( info != -7  )
    nbad++;

  fflush(stderr);
  mxsync_();
  if( me == 0 )
    fprintf( stderr, " \n pdspevx 8  invalid input \n");
  k = 3;
  ilb = -1;
  iub = 0;
  pdspevx(  &ivector, &k, &n, vecA, mapA, 
           &lb, &ub, &ilb, &iub, &abstol, &meigval, vecZ, mapZ, eval, iscratch,
          &isize, iptr, &ptr_size ,scratch, &rsize, &info);
  if( info != -8  )
    nbad++;

  fflush(stderr);
  mxsync_();
  if( me == 0 )
    fprintf( stderr, " \n pdspevx 9  invalid input \n");
  k = 3;
  ilb = 2;
  iub = 1;
  pdspevx(  &ivector, &k, &n, vecA, mapA, 
           &lb, &ub, &ilb, &iub, &abstol, &meigval, vecZ, mapZ, eval, iscratch,
          &isize, iptr, &ptr_size ,scratch, &rsize, &info);
  if( info != -9  )
    nbad++;

  fflush(stderr);
  mxsync_();
  if( me == 0 )
    fprintf( stderr, " \n pdspevx 9  invalid input \n");
  k = 3;
  ilb = 1;
  iub = n+1;
  pdspevx(  &ivector, &k, &n, vecA, mapA, 
           &lb, &ub, &ilb, &iub, &abstol, &meigval, vecZ, mapZ, eval, iscratch,
          &isize, iptr, &ptr_size ,scratch, &rsize, &info);
  if( info != -9  )
    nbad++;

  fflush(stderr);
  mxsync_();
  if( me == 0 )
    fprintf( stderr, " \n pdspevx 13 invalid input \n");
  mapZ[0] = -1;
  pdspevx(  &ivector, &irange, &n, vecA, mapA, 
           &lb, &ub, &ilb, &iub, &abstol, &meigval, vecZ, mapZ, eval, iscratch,
          &isize, iptr, &ptr_size ,scratch, &rsize, &info);
  if( info != -13 )
    nbad++;
  mapZ[0] = 0;

  fflush(stderr);
  mxsync_();
  if( me == 0 )
    fprintf( stderr, " \n pdspevx 16 invalid input \n");
  k = 1;
  pdspevx(  &ivector, &irange, &n, vecA, mapA, 
           &lb, &ub, &ilb, &iub, &abstol, &meigval, vecZ, mapZ, eval, iscratch,
          &k, iptr, &ptr_size ,scratch, &rsize, &info);
  if( info != -16 )
    nbad++;

  fflush(stderr);
  mxsync_();
  if( me == 0 )
    fprintf( stderr, " \n pdspevx 18 invalid input \n");
  fflush(stderr);
  k = 1;
  pdspevx(  &ivector, &irange, &n, vecA, mapA, 
           &lb, &ub, &ilb, &iub, &abstol, &meigval, vecZ, mapZ, eval, iscratch,
          &isize, iptr, &k,scratch, &rsize, &info);
  if( info != -18 )
    nbad++;

  fflush(stderr);
  mxsync_();
  if( me == 0 )
    fprintf( stderr, " \n pdspevx 20 invalid input \n");
  fflush(stderr);
  k = 1;
  pdspevx(  &ivector, &irange, &n, vecA, mapA, 
           &lb, &ub, &ilb, &iub, &abstol, &meigval, vecZ, mapZ, eval, iscratch,
          &isize, iptr, &ptr_size ,scratch, &k, &info);
  if( info != -20 )
    nbad++;

  if( nprocs > 1 ) {
    fflush(stderr);
    mxsync_();
    if( me == 0 )
      fprintf( stderr, " \n pdspevx 13 invalid input \n");
    fflush(stderr);
    
    for ( ii = 0;  ii < n; ii++ )
      mapZ[ii] = 0;
    pdspevx(  &ivector, &irange, &n, vecA, mapA, 
           &lb, &ub, &ilb, &iub, &abstol, &meigval, vecZ, mapZ, eval, iscratch,
            &isize, iptr, &ptr_size ,scratch, &rsize, &info);
  if( info != -13 )
    nbad++;
    for ( ii = 0;  ii < n; ii++ ) {
      k = ( ii % nprocs);
      mapZ[ii] = k;
    }
  }

/*--------------------------------------------------
 * Data differs on processors
 *--------------------------------------------------*/

  ndiff = nprocs / 2;

  if( me == ndiff && nprocs > 1 ) {

    fflush(stderr);
    mxsync_();
    fprintf( stderr, " \n pdspevx 1 on processor = %d  \n", ndiff);
    k = 0;
    pdspevx(  &k, &irange, &n, vecA, mapA,  
          &dk, &ub, &ilb, &iub, &abstol, &meigval, vecZ, mapZ, eval, iscratch,
            &isize, iptr, &ptr_size ,scratch, &rsize, &info);
    if( info != -51 )
      nbad++;

    fflush(stderr);
    fflush(stderr);
    mxsync_();
    fprintf( stderr, " \n pdspevx 6 on processor = %d  \n", ndiff);
    dk = -10.e0;
    pdspevx(  &ivector, &irange, &n, vecA, mapA,  
          &dk, &ub, &ilb, &iub, &abstol, &meigval, vecZ, mapZ, eval, iscratch,
            &isize, iptr, &ptr_size ,scratch, &rsize, &info);
    if( info != -51 )
      nbad++;

    fflush(stderr);
    mxsync_();
    fprintf( stderr, " \n pdspevx maps differ on processor = %d  \n", ndiff);
    mapA[0] = 1;
    mapZ[0] = 1;
    mapA[1] = 0;
    mapZ[1] = 0;
    pdspevx(  &ivector, &irange, &n, vecA, mapA, 
           &lb, &ub, &ilb, &iub, &abstol, &meigval, vecZ, mapZ, eval, iscratch,
            &isize, iptr, &ptr_size ,scratch, &rsize, &info);
  if( info != -51 )
    nbad++;
    mapA[0] = 0;
    mapZ[0] = 0;
    mapA[1] = 1;
    mapZ[1] = 1;

  } else if( nprocs > 1 ) {

    for ( i = 0; i < 3; i++ ) {
      fflush(stderr);
      mxsync_();
      pdspevx(  &ivector, &irange, &n, vecA, mapA, 
           &lb, &ub, &ilb, &iub, &abstol, &meigval, vecZ, mapZ, eval, iscratch,
              &isize, iptr, &ptr_size ,scratch, &rsize, &info);
  if( info != -51 )
    nbad++;
    }
  }
  return( nbad );
}
Integer tst_pdspev( n, vecA, mapA,  vecZ, mapZ, eval, iscratch,
                    isize, iptr, ptr_size ,scratch, rsize)

  Integer          n, *mapA, *mapZ, *iscratch,
                   isize, ptr_size, rsize;
  DoublePrecision  **vecA, **vecZ, **iptr;
  DoublePrecision  *eval, *scratch;

{
  Integer          k, ii, nbad, me, ndiff, nprocs, i, info;
  Integer          *inull;
  
  DoublePrecision  *dnull, **ddnull;

  extern Integer mxsync_(), mxmynd_(), mxnprc_();
  extern void    pdspev();

  me     = mxmynd_();
  nprocs = mxnprc_();

  nbad = 0;

/*--------------------------------------------------
 * Pass pointers to NULL.
 *--------------------------------------------------*/

  inull  = NULL;
  dnull  = NULL;
  ddnull = NULL;
  
  fflush(stderr);
  mxsync_();
  if( me == 0 )
    fprintf( stderr, " \n pdspev 1 = NULL \n");
  pdspev( inull, vecA, mapA, vecZ, mapZ, eval, iscratch,
          &isize, iptr, &ptr_size ,scratch, &rsize, &info);
  if( info != -1 )
    nbad++;

  fflush(stderr);
  mxsync_();
  if( me == 0 )
    fprintf( stderr, " \n pdspev 2 = NULL \n");
  pdspev( &n, ddnull, mapA, vecZ, mapZ, eval, iscratch,
          &isize, iptr, &ptr_size ,scratch, &rsize, &info);
  if( info != -2 )
    nbad++;

  fflush(stderr);
  mxsync_();
  if( me == 0 )
    fprintf( stderr, " \n pdspev 3 = NULL \n");
  pdspev( &n, vecA, inull, vecZ, mapZ, eval, iscratch,
          &isize, iptr, &ptr_size ,scratch, &rsize, &info);
  if( info != -3 )
    nbad++;

  fflush(stderr);
  mxsync_();
  if( me == 0 )
    fprintf( stderr, " \n pdspev 4  = NULL \n");
  pdspev( &n, vecA, mapA, ddnull, mapZ, eval, iscratch,
          &isize, iptr, &ptr_size ,scratch, &rsize, &info);
  if( info != -4  )
    nbad++;

  fflush(stderr);
  mxsync_();
  if( me == 0 )
    fprintf( stderr, " \n pdspev 5  = NULL \n");
  pdspev( &n, vecA, mapA, vecZ, inull, eval, iscratch,
          &isize, iptr, &ptr_size ,scratch, &rsize, &info);
  if( info != -5  )
    nbad++;

  fflush(stderr);
  mxsync_();
  if( me == 0 )
    fprintf( stderr, " \n pdspev 6  = NULL \n");
  pdspev( &n, vecA, mapA, vecZ, mapZ, dnull, iscratch,
          &isize, iptr, &ptr_size ,scratch, &rsize, &info);
  if( info != -6  )
    nbad++;

  fflush(stderr);
  mxsync_();
  if( me == 0 )
    fprintf( stderr, " \n pdspev 7  = NULL \n");
  pdspev( &n, vecA, mapA, vecZ, mapZ, eval, inull,
          &isize, iptr, &ptr_size ,scratch, &rsize, &info);
  if( info != -7  )
    nbad++;

  fflush(stderr);
  mxsync_();
  if( me == 0 )
    fprintf( stderr, " \n pdspev 8  = NULL \n");
  pdspev( &n, vecA, mapA, vecZ, mapZ, eval, iscratch,
          inull, iptr, &ptr_size ,scratch, &rsize, &info);
  if( info != -8  )
    nbad++;

  fflush(stderr);
  mxsync_();
  if( me == 0 )
    fprintf( stderr, " \n pdspev 9  = NULL \n");
  pdspev( &n, vecA, mapA, vecZ, mapZ, eval, iscratch,
          &isize, ddnull, &ptr_size ,scratch, &rsize, &info);
  if( info != -9  )
    nbad++;

  fflush(stderr);
  mxsync_();
  if( me == 0 )
    fprintf( stderr, " \n pdspev 10 = NULL \n");
  pdspev( &n, vecA, mapA, vecZ, mapZ, eval, iscratch,
          &isize, iptr, inull,scratch, &rsize, &info);
  if( info != -10 )
    nbad++;

  mxsync_();
  if( me == 0 )
    fprintf( stderr, " \n pdspev 11 = NULL \n");
  pdspev( &n, vecA, mapA, vecZ, mapZ, eval, iscratch,
          &isize, iptr, &ptr_size ,dnull, &rsize, &info);
  if( info != -11 )
    nbad++;

  fflush(stderr);
  mxsync_();
  if( me == 0 )
    fprintf( stderr, " \n pdspev 12 = NULL \n");
  pdspev( &n, vecA, mapA, vecZ, mapZ, eval, iscratch,
          &isize, iptr, &ptr_size ,scratch, inull, &info);
  if( info != -12 )
    nbad++;

  fflush(stderr);
  mxsync_();
  if( me == 0 )
    fprintf( stderr, " \n pdspev 13 = NULL \n");
  pdspev( &n, vecA, mapA, vecZ, mapZ, eval, iscratch,
          &isize, iptr, &ptr_size ,scratch, &rsize, inull );

  /* Can't check info since passed NULL for info. */

  fflush(stderr);
  mxsync_();
  if( me == 0 )
    fprintf( stderr, " \n pdspev 2 = NULL \n");
  dnull = vecA[0];
  vecA[0] = NULL;
  pdspev( &n, vecA, mapA, vecZ, mapZ, eval, iscratch,
          &isize, iptr, &ptr_size ,scratch, &rsize, &info );
  if( info != -2 )
    nbad++;
  vecA[0] = dnull;

  fflush(stderr);
  mxsync_();
  if( me == 0 )
    fprintf( stderr, " \n pdspev 4  = NULL \n");
  dnull = vecZ[0];
  vecZ[0] = NULL;
  pdspev( &n, vecA, mapA, vecZ, mapZ, eval, iscratch,
          &isize, iptr, &ptr_size ,scratch, &rsize, &info );
  if( info != -4  )
    nbad++;
  vecZ[0] = dnull;

/*--------------------------------------------------
 * Invalid input
 *--------------------------------------------------*/

  fflush(stderr);
  mxsync_();
  if( me == 0 )
    fprintf( stderr, " \n pdspev 1 invalid input \n");
  k = -1;
  pdspev( &k, vecA, mapA, vecZ, mapZ, eval, iscratch,
          &isize, iptr, &ptr_size ,scratch, &rsize, &info);
  if( info != -1 )
    nbad++;

  fflush(stderr);
  mxsync_();
  if( me == 0 )
    fprintf( stderr, " \n pdspev 3 invalid input \n");
  mapA[0] = -1;
  pdspev( &n, vecA, mapA, vecZ, mapZ, eval, iscratch,
          &isize, iptr, &ptr_size ,scratch, &rsize, &info);
  if( info != -3 )
    nbad++;
  mapA[0] = 0;

  fflush(stderr);
  mxsync_();
  if( me == 0 )
    fprintf( stderr, " \n pdspev 5  invalid input \n");
  mapZ[0] = -1;
  pdspev( &n, vecA, mapA, vecZ, mapZ, eval, iscratch,
          &isize, iptr, &ptr_size ,scratch, &rsize, &info);
  if( info != -5 )
    nbad++;
  mapZ[0] = 0;

  fflush(stderr);
  mxsync_();
  if( me == 0 )
    fprintf( stderr, " \n pdspev 8  invalid input \n");
  k = 1;
  pdspev( &n, vecA, mapA, vecZ, mapZ, eval, iscratch,
          &k, iptr, &ptr_size ,scratch, &rsize, &info);
  if( info != -8  )
    nbad++;

  fflush(stderr);
  mxsync_();
  if( me == 0 )
    fprintf( stderr, " \n pdspev 10 invalid input \n");
  fflush(stderr);
  k = 1;
  pdspev( &n, vecA, mapA, vecZ, mapZ, eval, iscratch,
          &isize, iptr, &k,scratch, &rsize, &info);
  if( info != -10 )
    nbad++;

  fflush(stderr);
  mxsync_();
  if( me == 0 )
    fprintf( stderr, " \n pdspev 12 invalid input \n");
  fflush(stderr);
  k = 1;
  pdspev( &n, vecA, mapA, vecZ, mapZ, eval, iscratch,
          &isize, iptr, &ptr_size ,scratch, &k, &info);
  if( info != -12 )
    nbad++;

  if( nprocs > 1 ) {
    fflush(stderr);
    mxsync_();
    if( me == 0 )
      fprintf( stderr, " \n pdspev 5  invalid input \n");
    fflush(stderr);
    
    for ( ii = 0;  ii < n; ii++ )
      mapZ[ii] = 0;
    pdspev( &n, vecA, mapA, vecZ, mapZ, eval, iscratch,
            &isize, iptr, &ptr_size ,scratch, &rsize, &info);
    if( info != -5  )
      nbad++;
    for ( ii = 0;  ii < n; ii++ ) {
      k = ( ii % nprocs);
      mapZ[ii] = k;
    }
  }

/*--------------------------------------------------
 * Data differs on processors
 *--------------------------------------------------*/

  ndiff = nprocs / 2;

  if( me == ndiff && nprocs > 1 ) {

    fflush(stderr);
    mxsync_();
    fprintf( stderr, " \n pdspev 1 on processor = %d  \n", ndiff);
    k = n+1;
    pdspev( &k, vecA, mapA, vecZ, mapZ, eval, iscratch,
            &isize, iptr, &ptr_size ,scratch, &rsize, &info);
    if( info != -51 )
      nbad++;

    fflush(stderr);
    mxsync_();
    fprintf( stderr, " \n pdspev maps differ on processor = %d  \n", ndiff);
    mapA[0] = 1;
    mapZ[0] = 1;
    mapA[1] = 0;
    mapZ[1] = 0;
    pdspev( &n, vecA, mapA, 
           vecZ, mapZ, eval, iscratch,
            &isize, iptr, &ptr_size ,scratch, &rsize, &info);
    if( info != -51 )
      nbad++;
    mapA[0] = 0;
    mapZ[0] = 0;
    mapA[1] = 1;
    mapZ[1] = 1;

  } else if( nprocs > 1 ) {

    for ( i = 0; i < 2; i++ ) {
      fflush(stderr);
      mxsync_();
      pdspev( &n, vecA, mapA, 
           vecZ, mapZ, eval, iscratch,
              &isize, iptr, &ptr_size ,scratch, &rsize, &info);
    if( info != -51 )
      nbad++;
    }
  }
  return( nbad );
}
