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
/* **************************************************************
 *
 *   C routine lsl_conjugation performs the conjugation of 
 * lower and 
 *
 *
 *
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "globalp.c.h"

#define max(a,b) ((a) > (b) ? (a) : (b))
#define min(a,b) ((a) < (b) ? (a) : (b))

void lsl_conjugation2 ( n, vecA, mapA, vecB, mapB, iwork, work, buff_ptr )
     Integer *n, *mapA, *mapB, *iwork;
     DoublePrecision **vecA, **vecB, *work, **buff_ptr;
     /*

       This routine performs the computation 
       
       A -> B.A.B^{T} 
       
       where A is a symmetric matrix stored in distributed column pack format
       using only its lower triangular form according to mapA, and B is a
       lower triangular matrix stored in distributed column pack storage
       format according to mapB
       
       The storage location of A is overwritten by B.A.B^{T}.
       The integer array mapA returns the location of the columns.
       
       Argument:
       
       n = dimension of the matrices ( assumed to be square )

       vecA[i] = memory location of the i-th vector that I own

       mapA[i] = processor id of processor which owns real column i

       mapB[i] = ditto as mapA

       vecB[i] = ditto as vecA

       iwork = integer of length p for broadcasting

       work = DoublePrecision scratch array of length 2*n
       ( assume that scratch is as large as n*(n+1)/2 -- this is a lot but ...)


       buffscratch = DoublePrecision precision storage; same size as that of A
       
       buff_ptr = array of pointers to DoublePrecision precision numbers; length same as mapvecA
       
       */
     
{
  Integer i, me, nvecsA, nvecsB, linfo, id;
  Integer *iscrat, *mapvecA, *mapvecB, *iptr;
  DoublePrecision *dscrat, **bufptr;
  extern Integer mxmynd_();
  extern Integer fil_mapvec_ ();

  extern void pmmlsl2();
  extern void lu_mxm();

/*
	extern void pmmLUL();
*/

  me = mxmynd_ ();

#ifdef DEBUG1
  fprintf(stderr, " in conjug22 me = %d \n", me);
#endif
  
  i = 0;
  me = mxmynd_ ();

  linfo = 0;
  iscrat = iwork;
  *iscrat = *n;
  iptr = iscrat + 1;

  iscrat = iwork;
  dscrat= work;
  bufptr = buff_ptr;

  mapvecA = iscrat;
  nvecsA = fil_mapvec_ ( &me, n, mapA, mapvecA );
  iscrat += nvecsA;
  mapvecB = iscrat;
  nvecsB = fil_mapvec_ ( &me, n, mapB, mapvecB );
  iscrat += nvecsB;
  
  id = 0;
  for ( i = 0; i < *n; i++ ){
    if ( mapA[i] != mapB[i] ) {
      id = 1;
      break;
    }
  }
  
/*
 *
 *   NOTE:  If the method of choosing pmmLSL,pmmLUL vs. pmmlsl2,lu_mxm2 is changed
 *          then this change must also be made in memreq.
 *
 */

  if ( id == 1 ) { 
    pmmLSL ( n, mapA, mapvecA, vecA, mapB, mapvecB, vecB, iscrat, dscrat, bufptr);
    pmmLUL ( n, mapA, mapvecA, vecA, mapB, mapvecB, vecB, iscrat, dscrat, bufptr);
  }
  else {
    
    /*
      de-symmetrize the matrix;
      all to all lower triangular multiplies;
      all to all upper triangular multiplies;
      */
    
#ifdef DEBUG
    for ( j = 0; j < nvecsB; j++ )
      fprintf(stderr, " me = %d mapvecB[ %d ] = %d \n", me, j, mapvecB[j]);
#endif

#ifdef DEBUG1
    fprintf(stderr, " just before pmmlsl2 \n");
#endif
        

    
    pmmlsl2 ( n, mapA, mapvecA, vecA, mapB, mapvecB, vecB, iscrat, dscrat, bufptr);

#ifdef DEBUG1
    fprintf(stderr, " exit pmmlsl2 \n");
#endif

    
#ifdef DEBUG
    for ( j = 0; j < nvecsB; j++ )
      fprintf(stderr, " me = %d mapvecB[ %d ] = %d \n", me, j, mapvecB[j]);
#endif

#ifdef DEBUG1
    fprintf(stderr, " just before lu_mxm2 \n");
#endif
    
    lu_mxm2( n, vecA, mapA, n, vecB, mapB, iscrat, dscrat);
    
#ifdef DEBUG1
    fprintf(stderr, " exit lu_mxm2 \n");
#endif

  }
  
#ifdef DEBUG
  fprintf(stderr, " out of  conjug22 me = %d \n", me);

  for ( j = 0; j < nvecsB; j++ )
    fprintf(stderr, " me = %d mapvecB[ %d ] = %d \n", me, j, mapvecB[j]);
  
  for ( i = 0; i < nvecsB; i++ )
    for( j = 0; j < *n-mapvecB[i]; j++ )
      fprintf(stderr, " me = %d vecB[ %d ][ %d ] = %g \n", me, i, j, vecB[i][j]);
  
  
#endif
  
  return;
}
