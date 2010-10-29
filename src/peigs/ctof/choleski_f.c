/*
 $Id$
 *
 *
 *
 *   This is the fortran code wrapper for Choleski factorization
 *  of a positive definite symmetric matrix A in compact column
 *  storage format.
 *
 *
 *      subroutine choleski ( n, colQ, mapQ, iwork, work, info)
 *                    integer n, mapQ(*), iwork(*), info
 *                    DoublePrecision precision colQ(*), work(*)
 *
 * Arguments:
 *
 *        n  = (input/integer) dimension of the matrix
 *
 *         colQ(*)  = (input/output/DoublePrecision precison array)
 *                    input matrix ( as described in the data structure
 *                    of the manual).
 *                    The output is the Choleski factor of the matrix
 *                    column distributed according to mapQ(*).
 *
 *         mapQ(*)  = (input/integer array)  length n, mapQ(i) = process
 *                             id of the processor that holds this column.
 *
 *         iwork(*) = (scratch/integer) length 2n
 *
 *         work(*)  = (scratch/DoublePrecision precision) length  2n + 1
 *                    ( actually nprocsQ + n + 1).
 *
 *                    Work must contain at least bufsiz bytes (see cmbbrf.h)
 *
 *         info     = (output/integer) = 0 ; exited normally
 *
 */

#include <stdio.h>

#include "globalp.c.h"

void choleski_( n, colQ, mapQ, iwork, work, info)
     Integer *n, *mapQ, *iwork, *info;
     DoublePrecision *colQ, *work;
{
  /*
    this is the Fortran 77 wrapper for the choleski factorization routine
    the matrix Q and W are assume to be in a column wrapped format

    */
  
  Integer icont,ibuz;
  Integer me, nvecsQ;
  Integer *iscrat, i, j, linfo;
  
  DoublePrecision **matQ;
  DoublePrecision *dscrat;
  
  extern Integer mxmynd_();
  extern Integer fil_mapvec_();
  extern void choleski();
  extern Integer count_list();
  extern void xerbla_(), pxerbla2_();
  extern void g_exit_();
  
   i = 0;
  me = mxmynd_();

  if( info == NULL ) {
    i = -6;
    xerbla_( "Choleski ", &i);
  }
  
  if( n == NULL ) {
    i = -1;
    xerbla_( "Choleski ", &i);
  }
  else
    if( *n < 1 ) {
      i = -1;
      xerbla_( "Choleski ", &i);
    }
    else
      {
	iscrat = mapQ;
	for( j = 0; j < *n; j++ ){
	  if( iscrat == NULL ){
	    i = -3;
	    xerbla_( "Choleski \n", &i);
	  }
	  else
	    iscrat++;
	}
      }
  
  /*
    
    at this point inputs are minimally acceptable for a given processor
    check to see if  n and mapQ are the same set of processors
    
    */
  
  iwork[0] = *n;
  iscrat = iwork + 1;
  for( i = 0; i < *n; i++ )
    *(iscrat++) = mapQ[i];
  
  linfo = 0;
  i = *n+1;
  iscrat = iwork + *n + 1;
  
  i = i * sizeof(Integer);
  pxerbla2_( &i, (char *) iwork, mapQ, n, iscrat, &linfo );
  
  g_exit_( &linfo, "Choleski:Mapping inconsistancies.\n", mapQ, n, iwork, work );
  
  /*
    if  everything seems to agree
    */
  
  me = mxmynd_();
  matQ =(DoublePrecision **) work;
  nvecsQ = count_list( me, mapQ, n);
  
  if( nvecsQ < 1) {
    *info = -50;
    return;
  }
  
  dscrat =(DoublePrecision *)( work + nvecsQ );
  icont=0;
     ibuz=0;
  for( i = 0; i < *n; i++ )
  {
    if( mapQ[i] == me){
    matQ[icont] = &colQ[ibuz];
     icont++;
      ibuz+= *n-i;
  }
   }
  
  iscrat  = iwork;
  choleski( n, matQ, mapQ, iscrat, dscrat, info );
  
  return;
}
