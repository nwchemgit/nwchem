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
#include <math.h>

#include "globalp.c.h"

#define max(a,b) ((a) > (b) ? (a) : (b))

void sonenrm ( n, colA, mapA, norm, iwork, work, info)
     Integer *n, *mapA, *iwork, *info;
     DoublePrecision **colA, *work, *norm;
{

/******************************************************
 *
 * C subroutine sonenrm
 *
 * This routines computes the one norm of a column
 * distributed symmetric matrix in packed storage format
 *
   Arguments
   ---------
   In the following:

     INTEGER          = "pointer to Integer"
     DOUBLE PRECISION = "pointer to DoublePrecision"

     me     = this processor's id (= mxmynd_())
     nprocs = number of allocated processors ( = mxnprc_())
     nvecsA = number of entries in mapA equal to me
                  (= count_list( me, mapA, n ))
     sDP    = sizeof( DoublePrecision )

       
   n....... (input) INTEGER
            size of the matrix A

   colA ... (input) array of pointers to DoublePrecision,
                    length (nvecsA)
            The part of matrix A owned by this processer stored
            in packed format, i.e., colA[i] points to the diagonal
            element of the i-th column (or equivalently row) of A
            owned by this processor, i = 0 to nvecsA-1.
                
   mapA ... (input) INTEGER array, length (n)
            The i-th column (or equivalently row) of A is 
            owned by processor mapA[i], i = 0 to n-1.

   norm ... (output) DOUBLE PRECISION
            The one-norm of A

   iwork... (workspace) INTEGER array, length( n+nvecsA )

   work.... (workspace) DOUBLE PRECISION array,
                        length( n + 1 + mxlbuf_() / sDP + 1 )
       
   info.... (output) INTEGER
            = 0, not currently used
 */
  
  static Integer IONE = 1;
  
  Integer ll, nprocs, i, me, nvecsA, *mapvecA;
  Integer *proclist, jj, k, ii;
  Integer *iscrat;
  
  DoublePrecision scl;
  DoublePrecision *normvec, *workMX;
  
  extern DoublePrecision dasum_ ();


  extern void gsum00();
  extern void fil_dbl_list ();
  extern Integer mxmynd_();
  extern Integer fil_mapvec_();
  extern Integer reduce_list2();
  extern void fil_dbl_lst();

  me = mxmynd_ ();

  *info = 0;
  ll = *n;
  *norm = 0.e0;

  iscrat   = iwork;

  mapvecA  = iscrat;
  nvecsA   = fil_mapvec_ ( &me, &ll, mapA, mapvecA);
  iscrat  += nvecsA;
  
  if ( nvecsA == 0 )
    return;
  
  proclist = iscrat;
  nprocs   = reduce_list2( *n, mapA, proclist);
  iscrat  += nprocs;
  
  normvec = work;
  workMX  = work + *n + 1;
  fil_dbl_lst ( *n, normvec, 0.0e0);    /* zero out normvec */
  
  for ( i = 0; i < nvecsA; i++ ) {
    jj = mapvecA[i];
    ii = ll - jj;
    scl = dasum_ ( &ii, colA[i], &IONE );
    normvec[jj] = scl;
  }
  
  for ( i = 0; i < nvecsA; i++ ) {
    jj = mapvecA[i];
    for ( k = 1; k < *n-jj; k++ )
       normvec[jj+k] += fabs( colA[i][k] );
  }
  
  gsum00( (char *) normvec, ll * sizeof(DoublePrecision), 5, 10, proclist[0], nprocs, proclist, workMX);
  
  scl = 0.0e0;
  for ( i = 0; i < ll; i++ )
    scl = max( normvec[i], scl);
  
  *norm = scl;

  return;
}
