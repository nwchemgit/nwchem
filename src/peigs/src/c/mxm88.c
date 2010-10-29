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

#include "globalp.c.h"

#define ZERO ((DoublePrecision) 0.0e0)

#define min(a,b) ((a) < (b) ? (a) : (b))

void mxm88 ( n, colQ, mapQ, m, colW, mapW, iwork, work, iptr)
     Integer *n, *mapQ, *m, *mapW, *iwork;
     DoublePrecision **colQ, **colW, *work, **iptr;

/**************************************************************
 *
 * Subroutine mxm88
 *
   This subroutine does the matrix multiplication W <- Q * W
       
   where Q is a n x n symmetric matrix in packed storage format,
   and row (or equivalently column) distributed,
   and matrix W is a general n x m matrix distributed by columns.
       
   ARGUMENTS
   ---------
   In the following:

     INTEGER          = "pointer to Integer"
     DOUBLE PRECISION = "pointer to DoublePrecision"

     me     = this processor's id (= mxmynd_())
     nprocs = number of allocated processors ( = mxnprc_())
     nvecsW = number of entries in mapW equal to me
                  (= count_list( me, mapW, m ))
     nvecsQ = number of entries in mapQ equal to me
                  (= count_list( me, mapQ, n ))
     sDP    = sizeof( DoublePrecision )

   n ...... (input) INTEGER
            n is the dimension of the symmetric matrix Q

   colQ ... (input) array of pointers to DoublePrecision,
                    length (nvecsQ)
            The part of matrix Q owned by this processer stored
            in packed format, i.e., colQ[i] points to the diagonal
            element of the i-th column (or equivalently row) of Q
            owned by this processor, i = 0 to nvecsQ-1.
                
   mapQ ... (input) INTEGER array, length (n)
            The i-th column (or equivalently row) of Q is 
            owned by processor mapQ[i], i = 0 to n-1.

   m ...... (input) INTEGER
            m is the number of columns in W.
       
   colW ... (input/output) array of pointers to DoublePrecision,
                           length (nvecsW)

            On Entry:
              The part of matrix W owned by this processer stored
              in packed format by columns, i.e., colW[i] points
              to the first element of the i-th column of W
              owned by this processor, i = 0 to nvecsW-1.
                
            On Exit:
              The result matrix Q * W stored in the same manner as
              W on entry.

   mapW ... (input) INTEGER array, length (m)
            The i-th column of W is 
            owned by processor mapW[i], i = 0 to m-1.
                

   iwork .. (workspace) INTEGER array, length ( 2*(n+m)+nvecsW+nvecsQ )

   work ... (workspace) DOUBLE PRECISION array, 
                        length ( maximum ( (nvecsW+1)*n, mxlbuf_()/sDP + 1 )

   iptr ... (workspace ) array of pointers to DoublePrecision,
                         length (nvecsW)

*/

{
  
  static Integer IONE = 1, MSGTYP = 51;
  
  Integer ll, k, i, *iscrat;
  Integer isize, me;
  Integer nvecsW, nvecsQ, j, iQ;
  Integer *mapvecQ, *mapvecW, *proclist, nprocs, npro, nele;
  Integer linfo;
  
  DoublePrecision t;
  DoublePrecision *buffer, scl;
  DoublePrecision *ptr, *ptr1;
  DoublePrecision  **Wptr;
  
  /*
    blas call
    */
  
  extern DoublePrecision ddot_();
  extern void pxerbla2_();
  extern void daxpy_();
  extern void dcopy_();
  extern void zero_out();
  
  extern void chol_pipe_bcast();
  extern Integer fil_mapvec_();
  extern Integer count_list();
  extern Integer reduce_list2();
  
  extern Integer mxwrit_(), mxread_();
  extern Integer menode_(), mxmynd_(), mxnprc_();
  
  extern void xerbla_();
  extern void g_exit_();
  extern void bbcast00();
  
  i = 0;
  me = mxmynd_();

  if ( m == NULL ) {
    i = -4;
    xerbla_( "MXM8 ", &i);
  }
  
  if ( n == NULL ) {
    i = -1;
    xerbla_( "MXM8 ", &i);
  }
  
  
  if ( *n < 1 ) {
    i = -1;
    xerbla_( "MXM8 ", &i);
  }
  
  iscrat = mapQ;
  for ( j = 0; j < *n; j++ ) {
    if ( iscrat++ == NULL ) {
      i = -3;
      xerbla_( "MXM8 \n", &i);
    }
  }
  iscrat = NULL;
  
  iscrat = mapW;
  for ( j = 0; j < *m; j++ ) {
    if ( iscrat++ == NULL ) {
      i = -6;
      xerbla_( "MXM8 \n", &i);
    }
  }
  iscrat = NULL;
  
  /*
    at this point inputs are minimally acceptable
    
    check to see if mapQ and mapW are the same set of processors
    */
  
  iscrat = iwork;
  j = 0;
  
/*
 * mapdif1_( n, mapQ, m, mapW, iscrat, &j );
 * 
 * if ( j != 0 ) {
 *   i = -3;
 * }
 * 
 */
  me = mxmynd_();
  nvecsQ = count_list( me, mapQ, n);
  nvecsW = count_list( me, mapW, m);
  
  for ( j = 0; j < nvecsQ; j++ )
    if ( colQ[j] == NULL ) {
      linfo = -2;
      i = min(linfo, i);
      break;
    }
  
  for ( j = 0; j < nvecsW; j++ )
    if ( colW[j] == NULL ) {
      i = min(i, -6);
      break;
    }
  
  g_exit_( &i, "Mapping problem or memory assignment problem in MXM8 \n", mapQ, n, iwork, work );
  
  linfo = 0;
  j = *n * sizeof(Integer);
  pxerbla2_( &j , (char *) mapQ, mapQ, n, iscrat, &i );
  linfo = min(linfo, i);
  j = *m * sizeof(Integer);;
  pxerbla2_( &j, (char *) mapW, mapQ, n, iscrat, &i);
  
  linfo = min(linfo, i);
  g_exit_( &linfo, "Mapping inconsistancies calling MXM8 \n", mapQ, n, iwork, work );
  
  ll = *n;
  
  /*
    should perform global sync to check for errors
    */
  
  me = mxmynd_();
/*
 * npro = mxnprc_();
 */
  npro = *m + *n;

  iscrat = iwork;

  proclist = iscrat;
  for ( i = 0; i < *m; i++ )
     iscrat[npro+i] = mapW[i];

  for ( i = 0; i < *n; i++ )
     iscrat[npro + *m + i] = mapQ[i];

  nele = *n + *m;
  nprocs = reduce_list2( nele, &iscrat[npro], proclist );
  iscrat += nprocs;
  
  mapvecQ = iscrat;
  nvecsQ = fil_mapvec_( &me, &ll, mapQ, mapvecQ );
  iscrat += nvecsQ;
  
  mapvecW = iscrat;
  nvecsW = fil_mapvec_( &me, m, mapW, mapvecW);
  iscrat += nvecsW;
  
  if (( nvecsQ == 0 ) && ( nvecsW == 0 ))
    return; /* no work */
  
  k = 0;
  Wptr = iptr;
  ptr = work;
  for ( i = 0; i < nvecsW; i++ ) {
    Wptr[i] = ptr;
    dcopy_( &ll, colW[i], &IONE, ptr, &IONE );
    ptr += ll;
  }
  
  buffer = ptr;
  
  /*
    zero out W
    */
  
  for ( i = 0; i < nvecsW; i++ ) {
    zero_out ( ll, colW[i] );
  }
  
  iQ = 0;
  for ( i = 0; i < ll; i++ ) {
    /*
      copy the column for passing
      */
    if ( mapQ[i] == me ) {
      isize = ll - i;
      dcopy_( &isize, colQ[iQ], &IONE, buffer, &IONE );
      
      /*
	broadcast to all processors holding W columns
	*/
      
      isize *= sizeof(DoublePrecision);
/*
 *    chol_pipe_bcast( buffer, isize, MSGTYP + i, mapQ[i], nprocs, proclist, iscrat);
 */
      bbcast00((char *) buffer, isize, MSGTYP + i, mapQ[i], nprocs, proclist);

      iQ++;
    }
    else {
      isize = (ll - i )* sizeof(DoublePrecision);
/*
 *     chol_pipe_bcast( buffer, isize, MSGTYP + i , mapQ[i], nprocs, proclist, iscrat);
 */
      bbcast00((char *) buffer, isize, MSGTYP + i, mapQ[i], nprocs, proclist);
    }
    
    /*
      the buffer contains the column Q[i:n-1, i]
      */
    
    for ( k = 0; k < nvecsW; k++ ) {
      isize = ll - i;
      ptr = Wptr[k] + i;
      t = *ptr;
      scl = ddot_( &isize, ptr, &IONE, buffer, &IONE);
      ptr1 = colW[k] + i;
      *ptr1 += scl;
      isize += (-1);
      daxpy_( &isize, &t, buffer+1 , &IONE, ptr1+1, &IONE);
    }
  }
  
  return;
}
