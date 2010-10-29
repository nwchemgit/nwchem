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
#include <memory.h>
#include <math.h>

#include "globalp.c.h"

#define max(a,b) ((a) > (b) ? (a) : (b))
#define min(a,b) ((a) < (b) ? (a) : (b))

/*
  Internal version of mxm5x which does not check mapdif.
  not supported for external use.
  
  this subroutine multiplies an upper triangular matrix with a full
  matrix the resulting matrix is stored in the full matrix 
  
  we assume that both L is an n by n matrix
  and U is an n by m matrix
  
  L is a lower triangular matrix;

  L <- L.U

  The product is stored in L.


  MUST have n = m

  Usage:

  Application:  L <- L.U
  
  */

void lu_mxm2( n, Lmatrix, mapL, m, colU, mapU, iscratch, scratch)
     Integer *n, *mapL, *mapU,  *m, *iscratch;
     DoublePrecision **colU, **Lmatrix, *scratch;
{
  /*
    n            = dimension of the matrix
    mapL         = index array holding the proces
    mapvecL      = index array holding the i-th column this processor owns
    colU         = DoublePrecision pointer to the array location of the i-th column
    
    ditto for mapvecL, Lmatrix, mapL
    scratch      = (DoublePrecision ) of size 2n, scratch buffer for message passing
                   scarch must contain at least bufsiz bytes (see cmbbrf.h)

    iscratch     = integer scratch space of size 2n
    
    remark: one can improve the x-fer by taking advantage size of the vector
    */
  
  static Integer IONE = 1;
  Integer jndx;
  Integer i, k, me, isize, *iptr, indx;
  Integer nvecsL, nvecsU;
  Integer j, linfo, ll, mm;
  Integer osize, rsize;		/* */
  Integer nproc, me_indx, last_proc , next_proc ;
  Integer *mapvecU, *mapvecL;
  Integer *iscrat, *proclist, i_dummy, iii;
  Integer *mapvec_in, maxsz;

  DoublePrecision *buffer, *dptr, *in_buffer, *out_buffer;
  DoublePrecision *Uvec, *Lvec, t;
  
  /*
    blas calls
    */

  extern void dcopy_(), daxpy_();
  extern DoublePrecision ddot_();
  extern Integer mxmynd_();
  extern void xerbla_();
  extern void mapdif1_();
  extern void gshellsort_();
  extern Integer indxL();

  /* 
    mxsubs calls
    */

  extern Integer mxwrit_(),  mxread_();
  extern Integer count_list();
  extern Integer fil_mapvec_();
  extern Integer ci_size_();
  extern void zero_out();
  extern Integer reduce_list2();
  
  i = 0;
  me = mxmynd_();

#ifdef DEBUG1
  fprintf(stderr, "node %d  In lu_mxm2 ****************************_____________ \n", me );
#endif

  if ( m == NULL )
    i = -4;
  else
    if ( n == NULL ) {
      i = -1;
      xerbla_( "lu_mxm2 ", &i);
    }
    else if ( *n < 1 )
      i = -1;
    else if ( *n != *m  )
      i = -3;
    else {
      iscrat = mapL;
      for ( j = 0; j < *n; j++ ) {
	if ( iscrat == NULL ) {
	  i = -3;
	  xerbla_( "lu_mxm2 \n", &i);
	}
	else
	  iscrat++;
      }
      iscrat = mapU;
      for ( j = 0; j < *m; j++ ) {
	if ( iscrat == NULL ) {
	  i = -6;
	  xerbla_( "lu_mxm2 \n", &i);
	}
	else
	  iscrat++;
      }
    }

  /*
    at this point inputs are minimally acceptable

    check to see if mapL and mapF are the same set of processors
    */

  iscrat = iscratch;
  mapdif1_( n, mapL, m, mapU, iscrat, &j );

  /*
   * if ( j != 0 ) {
   *   i = -3;
   * }
   */
  
  
  me = mxmynd_();
  nvecsL = count_list ( me, mapL, n);
  nvecsU = count_list ( me, mapU, m);
  
  if ( nvecsL + nvecsU == 0 )
    return;
  
#ifdef DEBUG
  fprintf (stderr, "me = %d, i = %d, n=%d, m=%d \n", me, i, *n, *m);
  for ( j = 0; j < *n; j++ )
    fprintf (stderr, "lu_mxm21 me = %d, mapL[%d]= %d, mapU[%d]=%d, \n", me, j, mapL[j], j, mapU[j]);
#endif

  
  for ( j = 0; j < nvecsL; j++ )
    if ( Lmatrix[j] == NULL ) {
      linfo = -2;
      i = min(linfo, i);
      break;
    }

  for ( j = 0; j < nvecsU; j++ )
    if ( colU[j] == NULL ) {
      i = min(i, -6);
      break;
    }

  /*

    g_exit_( &i, "Mapping problem or memory assignment problem lu_mxm2X \n", mapL, n, iscratch, scratch );

    linfo = 0;
    ll = *n * sizeof(Integer);
    pxerbla2_( &ll, mapL, mapL, n, iscrat, &i );
    linfo = min(linfo, i);
    ll = *m * sizeof(Integer);

    pxerbla2_( &ll, mapU, mapL, n, iscrat, &i);

    linfo = min(linfo, i);

    g_exit_( &linfo, "Mapping inconsistancies calling lu_mxm2X \n", mapL, n, iscratch, scratch );

    */

  me = mxmynd_();
  
  ll = *n;
  mm = *m;

#ifdef DEBUG1
  fprintf(stderr, " 1 in lu_mxm2.c me = %d \n", me);
#endif
  
  
  iscrat = iscratch;
  mapvecL = iscrat;
  nvecsL = fil_mapvec_( &me, &ll, mapL, mapvecL );
  iscrat += nvecsL;
  mapvecU = iscrat;
  nvecsU = fil_mapvec_( &me, &mm, mapU, mapvecU );
  iscrat += nvecsU;
  
#ifdef DEBUG1
  fprintf(stderr, " 2 in lu_mxm2.c me = %d \n", me);
#endif
  
  iptr = (Integer *) scratch;
  for ( i = 0; i < ll; i++ )
    *(iptr++) = mapL[i];
  for ( i = 0; i < mm; i++ )
    *(iptr++) = mapU[i];
  
#ifdef DEBUG1
  fprintf(stderr, " 3 in lu_mxm2 1 me = %d \n", me);
#endif
  
  iptr = (Integer *) scratch;
  indx = ll + mm;
  
  nproc = reduce_list2( indx, iptr, iscrat);
  
  proclist = iscrat;
  iscrat += nproc;
  mapvec_in = iscrat;
  
  gshellsort_( &nproc, proclist);
  
  indx = indxL( me, nproc, proclist);
  me_indx = indx;
  indx += nproc;
  last_proc = (indx - 1) % nproc;
  last_proc = proclist[last_proc];
  next_proc = (indx + 1) % nproc;
  next_proc = proclist[next_proc];
  

#ifdef DEBUG1
  fprintf (stderr, "4 me = %d, last_proc = %d, next_proc =%d, \n", me, last_proc, next_proc);
#endif
  
  maxsz = 0;
  isize = 0;
  for ( indx = 0; indx < nproc; indx++ ) {
    isize = ci_size_( &proclist[indx], &ll, mapL);
    isize *= 2;
    maxsz = max( isize, maxsz );
  }
  
#ifdef DEBUG1
  fprintf (stderr, "5 after maxsz me = %d, i = %d, n=%d, m=%d \n", me, i, *n, *m);
#endif
  
  buffer = (DoublePrecision *) scratch;
  in_buffer = buffer;
  indx = 0;
  for ( jndx = 0; jndx < nvecsU; jndx++ ){
    isize = ll - mapvecU[jndx];
    dcopy_( &isize, colU[jndx], &IONE, &buffer[indx], &IONE);
    indx += isize;
    isize = ll - mapvecL[jndx];
    dcopy_( &isize, Lmatrix[jndx], &IONE, &buffer[indx], &IONE);
    indx += isize;
  }
  buffer += indx;
  
  /*
    zero out the matrix
    */
  
  
#ifdef DEBUG
  fprintf (stderr, "6 after maxsz me = %d, i = %d, n=%d, m=%d \n", me, i, *n, *m);
#endif
  
  for ( jndx = 0; jndx < nvecsL; jndx++ ){
    isize = ll - mapvecL[jndx];
    zero_out ( isize, Lmatrix[jndx]);
  }
  
#ifdef DEBUG
  fprintf (stderr, "7 after maxsz me = %d, i = %d, n=%d, m=%d \n", me, i, *n, *m);
#endif
  
  out_buffer = in_buffer + maxsz;
  
  /*
    compute local matrix multiply
    */

  dptr = in_buffer;
  for ( jndx = 0; jndx < nvecsU; jndx++ ) {
    j = mapvecU[jndx];
    Uvec = dptr;
    Lvec = dptr + ll - j;
    for ( indx = 0; indx < nvecsU; indx++ ) {
      k = mapvecL[indx];
      if ( j <= k ){
	isize = ll - k;
	iii = k - j;
	t = Uvec[iii];
	daxpy_( &isize, &t, Lvec+iii, &IONE, Lmatrix[indx], &IONE);
      }
    }
    dptr += 2*(ll-j);
  }
  
#ifdef DEBUG
  fprintf (stderr, "9 after maxsz me = %d, i = %d, n=%d, m=%d \n", me, i, *n, *m);
#endif
  
  /*
    buffer     = out_buffer + maxsz;
    */
  
  isize = 2*ci_size_( &me, &ll, mapL);   /* the number Q vectors that I own */
  osize = isize;
  me_indx += nproc;
  
  for ( i = 0; i < nproc-1 ; i++ ) {
    dcopy_(&osize, in_buffer, &IONE, out_buffer, &IONE);
    if ( (me_indx % 2) ==  0 ) {
      rsize = osize * sizeof(DoublePrecision);
      if ( rsize != 0 ) {
	rsize = mxwrit_( out_buffer, &rsize, &next_proc, &i );
      }
      
      indx = ( me_indx-i-1 ) % nproc;
      indx = proclist[indx]; /* the current processor id number */
      nvecsL = fil_mapvec_( &indx, &ll, mapL, mapvec_in);
      isize = 2*ci_size_( &indx, &ll, mapL);   /* the number Q vectors that I own */
      
      /*
	 mapvec_in contains the mapvecQ for processor proclist[i])
	 */
      
      rsize = isize*sizeof(DoublePrecision);
      if ( rsize != 0 ){
	rsize = mxread_( in_buffer, &rsize, &last_proc, &i );
      }
    }
    else {
      indx = ( me_indx-i-1 ) % nproc;
      indx = proclist[indx]; /* the current processor id number */
      nvecsL = fil_mapvec_( &indx, &ll, mapL, mapvec_in);
      isize = 2*ci_size_( &indx, &ll, mapL);   /* the number Q vectors that I own */
      
      /*
	mapvec_in contains the mapvecQ for processor proclist[i])
	*/
      
      rsize = isize * sizeof(DoublePrecision);
      if ( rsize != 0 )
	rsize = mxread_( in_buffer, &rsize, &last_proc, &i );
      
      rsize = osize * sizeof(DoublePrecision);
      if ( rsize != 0 ){
	rsize = mxwrit_( out_buffer, &rsize, &next_proc, &i );
      }
    }
    
    osize = isize;
    
    dptr = in_buffer;
    for ( jndx = 0; jndx < nvecsL; jndx++ ) {
      j = mapvec_in[jndx];
      Uvec = dptr;
      Lvec = dptr + ll - j;
      for ( indx = 0; indx < nvecsU; indx++ ) {
	k = mapvecL[indx];
	if ( j <= k ){
	  i_dummy = ll - k;
	  iii = k - j;
	  t = Uvec[iii];
	  daxpy_( &i_dummy, &t, Lvec+iii, &IONE, Lmatrix[indx], &IONE);
	}
      }
      dptr += 2*(ll-j);
    }
  }
  
  return;
}





