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
/******************************************************
 *
 *  C subroutine onenrm
 *
 *
 *    computes the one norm of a column distributed
 *    n x n  matrix
 *
 *    void onenrm ( n, colA, mapA, norm, iwork, work, info)
 *           Integer *n, *mapA, *iwork, *info;
 *           DoublePrecision **colA, *work, *norm;
 *
 *   
 *    output: one norm
 *
 *
 */

#include <stdio.h>
#include <math.h>

#include "globalp.c.h"

#define max(a,b) ((a) > (b) ? (a) : (b))

/******************************************************
 *
 *  C subroutine onenrm
 *
 *
 *    computes the one norm of a column distributed
 *    n x n  matrix
 *
 *    void onenrm ( n, colA, mapA, norm, iwork, work, info)
 *           Integer *n, *mapA, *iwork, *info;
 *           DoublePrecision **colA, *work, *norm;
 *
 *   
 *    output: one norm
 *
 *
 */

void one_nrm ( n, m, colA,  mapA, norm, iwork, work)
     Integer *n, *m, *mapA, *iwork;
     DoublePrecision **colA, *norm, *work;
{
  /*
    computes the infinity norm of matrix A
    
    n       = (input/integer) vector length
    m       = (input/integer) number of vectors
    
    colA    = (input/1-D array of pointers to doubles ) pointer to the columns
    
    mapA    = (input/1-D array of integers) mapA[i] = the processor
    id which owns column i
    
    iwork   = (scratch space/integer) length  nvecsA + nprocs
    
    work    = DoublePrecision precision scratch space of size
    
    */
  
  static Integer IONE = 1;
  Integer ll, nprocs, i, me, nvecsA, *mapvecA;
  Integer *proclist;
  Integer *iscrat;
  DoublePrecision t, anorm;
  
  extern DoublePrecision dasum_();

  extern void gmax00();
  extern void xerbla_();
  extern Integer mxmynd_();
  extern Integer reduce_list2();
  extern Integer fil_mapvec_();
  extern Integer reduce_list22();
  
  me       = mxmynd_();
  
#ifdef DEBUG1
  fprintf(stderr, " me = %d n = %d, m = %d \n", me, *n, *m );
#endif
  if ( n == NULL ) {
    ll = -1;
    xerbla_("ONENRM\n", &ll);
  }
  
  ll = *m;
  iscrat   = iwork;
  mapvecA  = iscrat;
  nvecsA   = fil_mapvec_( &me, &ll, mapA, mapvecA);
  iscrat  += nvecsA;
  
  if ( nvecsA <= 0 )
    return;
  
  proclist = iscrat;
  nprocs   = reduce_list2( ll, mapA, proclist);
  iscrat  += nprocs;
  
#ifdef DEBUG1
  fprintf(stderr, " me = %d nvecsA = %d, nprocs = %d \n", me, nvecsA, nprocs);
  fprintf(stderr, " me = %d n = %d, m = %d \n", me, *n, *m );
  for ( i = 0; i < *m; i++ )
     fprintf(stderr, " me = %d mapA[%d] = %d \n", me, i, mapA[i]);
  for ( i = 0; i < nprocs; i++ )
     fprintf(stderr, " me = %d proclist[%d] = %d \n", me, i, proclist[i]);
#endif

  anorm = 0.0e0;
  for ( i = 0; i < nvecsA; i++ ) {
    t = dasum_( n, colA[i], &IONE);
    anorm = max(anorm, t);
  }
  
  /*
    spanning tree 
    */
  
  gmax00( (char *) &anorm, 1, 5, 16, proclist[0], nprocs, proclist, work);
  
  /*
    global combine max
    */
  
  
  *norm = anorm;
  return;
}


