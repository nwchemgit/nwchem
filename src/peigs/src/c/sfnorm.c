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
 *  C subroutine sonenrm
 *
 *
 *    computes the Frobenious(F)-norm of a column distributed symmetric
 *    matrix in packed storage format
 *
 *    void sfnorm ( n, mapA, colA, norm, iwork, work, info)
 *           Integer *n, *mapA, *iwork, *info;
 *           DoublePrecision **colA, *norm, *work;
 *
 *   
 *    output: F-norm of a symmetrix matrix distributed using
 *            a packed storage format
 *
 *
 */

#include <stdio.h>
#include <math.h>

#include "globalp.c.h"

void sfnorm( n, colA, mapA, norm, iwork, work, info)
     Integer *n, *mapA, *iwork, *info;
     DoublePrecision **colA, *work, *norm;
{
  /*
    computes the F-norm of a symmetric n-by-n matrix A, where

    F-norm(A) = sqrt( sum_{i,j} a_{i,j}^2 )
    
    n       = size of the matrix
    colA    = DoublePrecision pointer to the columns
    mapA    = distribution of the columns
    iwork   = integer scratch space
    work    = DoublePrecision precision scratch space
    
    */
  
  Integer ll, nprocs, i, me, nvecsA, *mapvecA, kk;
  Integer *proclist, jj, ii;
  Integer *iscrat;
  
  DoublePrecision root2, normvec[1], dummy_vec;
  DoublePrecision s, t, y, f, c; 
  extern DoublePrecision ddot_();
  
  extern void gsum00();
  extern Integer mxmynd_();
  extern Integer fil_mapvec_();
  extern Integer reduce_list2();
  extern void fil_dbl_lst();

  me = mxmynd_ ();
  ll = *n;
  iscrat   = iwork;
  mapvecA  = iscrat;
  nvecsA   = fil_mapvec_ ( &me, &ll, mapA, mapvecA);
  iscrat += nvecsA;
  
  /*
    
    no vectors here; just signal return; if no in mapA there's no way to
    get info unless the user explicitely does it some other way
    
    */
  
  *info = 0;
  if ( nvecsA == 0 ) {
    *info = -50;
    return;
  }
  
  iscrat  += nvecsA;
  proclist = iscrat;
  nprocs   = reduce_list2( *n, mapA, proclist);
  iscrat  += nprocs;
  
  dummy_vec = 0.0e0;
  s = 0.0e0;
  c = 0.0e0;
/*
  should keep the sum = sum_{i}^{N} ( 1 + epsilon(i)) x(i)
  with abs(epsilon) < 5*u + O(Nu**2) where  t = # mantissa digits

  u = machine unit ~  .5*Base^(1-t) if rounding 
  u = machine unit ~     Base^(1-t) if chopping
*/
  for ( i = 0; i < nvecsA; i++ ) {
    root2 = colA[i][0];
    y = c + root2*root2;
    t = s + y;
    f = 0.e0;
    if (( y < 0.0 && s < 0.0 ) || ( y > 0.0 && s > 0.0 ))
      f = ( 0.46*t - t ) + t;
    c = ((s-f)-(t-f)) + y;
    s = t;
  }
  s += c;
  normvec[0] = dummy_vec = s;
  
  dummy_vec = (DoublePrecision) 0.0e0;
  s = 0.0e0;
  c = 0.0e0;
  for ( i = 0; i < nvecsA; i++ ) {
    jj = mapvecA[i];
    ii = ll - jj - 1;
    for ( kk = 1; kk <= ii; kk++ ) {
      y = c + colA[i][kk]*colA[i][kk];
      t = s + y;
      f = 0.e0;
      if (( y < 0.0 && s < 0.0 ) || ( y > 0.0 && s > 0.0 ))
	f = ( 0.46*t - t ) + t;
      c = ((s-f)-(t-f)) + y;
      s = t;
    }
    s += c;
  }
  
  /*
    dummy_vec += ddot_( &ii, &colA[i][1], &IONE, &colA[i][1], &IONE);
    */
  
  dummy_vec = s;
  dummy_vec *= (DoublePrecision) 2.0e0;
  normvec[0] += dummy_vec;
  
#ifdef DEBUG
  for ( ii = 0; ii < nprocs; ii++ )
    fprintf(stderr, " me = %d proclist[%d] = %d \n", me, ii, proclist[ii]);
#endif
  
  gsum00( (char *) normvec, 1, 5, 10, proclist[0], nprocs, proclist, work);
  
#ifdef DEBUG
  for ( ii = 0; ii < *n; ii++ )
    fprintf(stderr, " 2 me = %d normvec[%d] = %g \n", me, ii, normvec[ii]);
#endif
  
  normvec[0] = sqrt(normvec[0]);
  *norm = normvec[0];
  
#ifdef DEBUG
  fprintf(stderr, " sfnorm: me = %d norm = %g \n", me, *norm );
#endif
  
  return;
}


