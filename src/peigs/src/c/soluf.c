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

/*
  PeIGS internal routine, currently unsupported.
  
  Let U = an upper triangular matrix
  Let W = full matrix
  
  This subroutine solves the matrix Y for the matrix problem
  
  U Y = W
  
  The matrix Y overwrites W.
  
  Data structure:
  
  The matrix U is row wrapped.
  The matrix W is column wrapped.
  
  */

void upperUF_ ( n, rowU, mapU, nW, colW, mapW, iwork, work)
     Integer *n, *mapU, *nW, *mapW, *iwork;
     DoublePrecision **rowU, **colW, *work;
{
  /*
    n            = dimension of the matrix
    mapU         = index array holding the proces
    rowU         = DoublePrecision pointer to the array location of the i-th column
    
    ditto for mapvecW, colW, mapW
    
    iwork = integer scratch space of size p
    
    remark: one can improve the x-fer by taking advantage size of the
    vector by segmenting the column vectors
    
    should put an info in the argument list with regard to underflow
    
    */
  
  static Integer IONE = 1;
  
  Integer i, iii, k, me, isize;
  Integer nvecsU, nvecsW;
  Integer count_list(), i_U;
  Integer nprocs;
  Integer *mapvecU, *mapvecW;
  Integer *iscrat;
  Integer *proclist;

  DoublePrecision *ptr;
  DoublePrecision t, *buffer;
  DoublePrecision dummy;
  DoublePrecision safeulp;
  
  extern void combine_vector();
  extern Integer fil_mapvec_();
  extern Integer indxL();
  
  /*
    blas calls
    */
  
  extern void dcopy_ ();
  extern DoublePrecision ddot_ ();

  
  /*
    mxsubs communication calls
    */
  
  extern Integer mxwrit_ (), mxread_ (), mxmynd_ ();
  extern void chol_pipe_bcast();
  extern Integer mxmynd_ (), reduce_list2();
  
  me = mxmynd_ ();
  
  /*
    info = 0;
   */

  safeulp = DLAMCHS;

  iscrat = iwork;
  mapvecU = iscrat;
  nvecsU = fil_mapvec_ ( &me, n, mapU, mapvecU );
  
  iscrat += nvecsU;
  mapvecW = iscrat;
  nvecsW = fil_mapvec_ ( &me, nW, mapW, mapvecW );
  iscrat += nvecsW;

  proclist = iscrat;
  nprocs = reduce_list2( *nW, mapW, proclist);
  iscrat += nprocs;
  
  if ( nvecsU + nvecsW == 0 )
    return;
  
  buffer = work;
  
  /*
    count number of columns that I own
    */
  
  
  i_U = nvecsU-1;
  for ( i = *n - 1 ; i >= 0; i-- ) {
    isize = *n - i;
    if ( mapU[i] == me ) {
      /*
	indx = indxL ( i, nvecsU, mapvecU );
	*/
      dcopy_ ( &isize, rowU[i_U], &IONE, buffer, &IONE );
      isize *= sizeof(DoublePrecision);
      chol_pipe_bcast( (char *) buffer, isize, i, mapU[i], nprocs, proclist, iscrat );
      i_U--;
    }
    else {
      isize = *n - i;
      isize *= sizeof(DoublePrecision);
      iii = mapU[i];
      chol_pipe_bcast( (char *) buffer, isize, i, iii, nprocs, proclist, iscrat);
    }
    
    /*
      the buffer[i:n-1] now contains U[i,i:n-1]
      */

   *buffer = (DoublePrecision) 1.0 / *buffer;

/*
    if ( *buffer > safeulp )
      *buffer = (DoublePrecision) 1.0 / *buffer;
      else {
      *info = 1;
      *buffer = 1.0e0/safeulp;
    }
    
*/
    dummy = *buffer;
    isize = *n-i-1;
    
    for ( k = 0; k < nvecsW; k++ ) {
      ptr = &colW[k][i + 1];
      t = (DoublePrecision) -ddot_ ( &isize, &buffer[1], &IONE, ptr, &IONE);
      t += colW[k][i];
      t *= dummy;
      colW[k][i] = t;
    }
  }
  
  return;
}


      
