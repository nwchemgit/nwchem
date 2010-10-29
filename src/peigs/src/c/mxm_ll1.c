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
  and F is an n by m matrix
  
  L is a lower triangular matrix;

  F <- L F

  The lower triangular matrix L is stored in distributed column format.
  The full matrix F is stored in distributed column format.  No assumption
  is made about its structure.

  The product is stored in F.

  Usage:

  Application:  F <- L.F

  */

void mxm_llx( n, colL, mapL, m, colF, mapF, iscratch, scratch)
     Integer *n, *mapL, *mapF,  *m, *iscratch;
     DoublePrecision **colF, **colL, *scratch;
{
  /*
    n            = dimension of the matrix
    mapL         = index array holding the proces
    mapvecL      = index array holding the i-th column this processor owns
    colU         = DoublePrecision pointer to the array location of the i-th column
    
    ditto for mapvecL, colL, mapL
    scratch      = (DoublePrecision ) of size 2n, scratch buffer for message passing
                   scratch must contain at least bufsiz bytes (see cmbbrf.h)

    iscratch     = integer scratch space of size 2n
    
    remark: one can improve the x-fer by taking advantage size of the vector
    */
  
  static Integer IONE = 1;
  Integer jndx;
  Integer i, k, me, isize, *iptr, indx;
  Integer nvecsL, nvecsF, ii;
  Integer j, linfo, ll, mm, fff ;
  Integer nproc;
  Integer osize, rsize;		/* */

  Integer *mapvecF, *mapvecL;
  Integer *iscrat, *proclist;
  Integer *mapvec_in, me_indx, last_proc, next_proc, maxsz;

  DoublePrecision *buffer, *d_ptr, *in_buffer, *out_buffer;
  DoublePrecision *F_ptr, *L_ptr, t;

  /*
    blas calls
    */

  extern void dcopy_ (), daxpy_ ();
  extern DoublePrecision ddot_ ();

  /* 
    mxsubs calls
    */

  extern Integer mxwrit_ (),  mxread_ ();
  extern Integer count_list();
  extern Integer fil_mapvec_ ();
  extern Integer ci_size_ ();

  i = 0;
  me = mxmynd_ ();

  if ( m == NULL )
    i = -4;
  else
    if ( n == NULL ) {
      i = -1;
      xerbla_ ( "mxm_ll ", &i);
    }
    else if ( *n < 1 )
      i = -1;
    else {
      iscrat = mapL;
      for ( j = 0; j < *n; j++ ) {
	if ( iscrat == NULL ) {
	  i = -3;
	  xerbla_ ( "mxm_ll \n", &i);
	}
	else
	  iscrat++;
      }
      iscrat = mapF;
      for ( j = 0; j < *m; j++ ) {
	if ( iscrat == NULL ) {
	  i = -6;
	  xerbla_ ( "mxm_ll \n", &i);
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
  mapdif1_ ( n, mapL, m, mapF, iscrat, &j );

  /*
   * if ( j != 0 ) {
   *   i = -3;
   * }
   */
  
  
  me = mxmynd_ ();
  nvecsL = count_list ( me, mapL, n);
  nvecsF = count_list ( me, mapF, m);
  
  if ( nvecsL + nvecsF == 0 )
    return;
  
#ifdef DEBUG
  fprintf (stderr, "me = %d, i = %d, n=%d, m=%d \n", me, i, *n, *m);
  for ( j = 0; j < *n; j++ )
    fprintf (stderr, "mxm_ll1 me = %d, mapL[%d]= %d, mapF[%d]=%d, \n", me, j, mapL[j], j, mapF[j]);
#endif

  
  for ( j = 0; j < nvecsL; j++ )
    if ( colL[j] == NULL ) {
      linfo = -2;
      i = min(linfo, i);
      break;
    }

  for ( j = 0; j < nvecsF; j++ )
    if ( colF[j] == NULL ) {
      i = min(i, -6);
      break;
    }

  /*

    g_exit_ ( &i, "Mapping problem or memory assignment problem mxm_llX \n", mapL, n, iscratch, scratch );
    g_exit_ ( &i, "Mapping problem or memory assignment problem mxm_ll \n", mapF, n, iscratch, scratch );

    linfo = 0;
    ll = *n * sizeof(Integer);
    pxerbla2_ ( &ll, mapL, mapL, n, iscrat, &i );
    linfo = min(linfo, i);
    ll = *m * sizeof(Integer);

    pxerbla2_ ( &ll, mapF, mapL, n, iscrat, &i);

    linfo = min(linfo, i);

    g_exit_ ( &linfo, "Mapping inconsistancies calling mxm_llX \n", mapL, n, iscratch, scratch );
    g_exit_ ( &linfo, "Mapping inconsistancies calling mxm_llX \n", mapF, n, iscratch, scratch );

    */

  me = mxmynd_ ();

  ll = *n;
  mm = *m;
  
  iscrat = iscratch;
  mapvecL = iscrat;
  
#ifdef DEBUG1
  fprintf(stderr, " 2 in mxm_llx me = %d \n", me);
#endif
  
  nvecsL = fil_mapvec_ ( &me, &ll, mapL, mapvecL );
  iscrat += nvecsL;

#ifdef DEBUG1
  fprintf(stderr, " 22 in mxm_llx me = %d \n", me);
#endif
  
  mapvecF = iscrat;
  nvecsF = fil_mapvec_ ( &me, m, mapF, mapvecF );
  iscrat += nvecsF;
  
  iptr = (Integer *) scratch;
  for ( i = 0; i < ll; i++ )
    *(iptr++) = mapL[i];
  for ( i = 0; i < mm; i++ )
    *(iptr++) = mapF[i];
  
#ifdef DEBUG1
  fprintf(stderr, " 3 in mxm_llx 1 me = %d \n", me);
#endif

  iptr = (Integer *) scratch;
  nproc = reduce_list2( ll+mm, iptr, iscrat);
  proclist = iscrat;
  iscrat += nproc;
  mapvec_in = iscrat;

  gshellsort_ ( &nproc, proclist);
  
  indx = indxL ( me, nproc, proclist);
  me_indx = indx;
  last_proc = (indx + nproc - 1) % nproc;
  last_proc = proclist[last_proc];
  next_proc = (indx + nproc + 1) % nproc;
  next_proc = proclist[next_proc];
  
#ifdef DEBUG1
  fprintf (stderr, "4 me = %d, i = %d, n=%d, m=%d \n", me, i, *n, *m);
#endif
  
  maxsz = 0;
  isize = 0;
  for ( indx = 0; indx < nproc; indx++ ) {
    isize = ci_size_ ( &proclist[indx], &ll, mapL);
    /*
      the number Q vectors that I own
      */
    maxsz = max( isize, maxsz );
  }
  
#ifdef DEBUG1
  fprintf (stderr, "5 after maxsz me = %d, i = %d, n=%d, maxsz=%d \n", me, i, *n, maxsz);
#endif
  
  buffer = (DoublePrecision *) scratch;
  indx = 0;
  for ( jndx = 0; jndx < nvecsF; jndx++ ){
    dcopy_ ( &ll, colF[jndx], &IONE, &buffer[indx], &IONE);
    indx += ll;
  }
  buffer += indx;
  
  /*
    zero out the matrix
    */
  
#ifdef DEBUG1
  fprintf (stderr, "6 after maxsz me = %d, i = %d, n=%d, m=%d \n", me, i, *n, *m);
#endif
  
  for ( jndx = 0; jndx < nvecsF; jndx++ ){
    d_ptr = colF[jndx];
    fil_dbl_lst ( ll, d_ptr, 0.0e0);
  }

#ifdef DEBUG1
  fprintf (stderr, "7 after maxsz me = %d, i = %d, n=%d, m=%d \n", me, i, *n, *m);
#endif


  in_buffer = buffer;
  out_buffer = buffer + maxsz;

  /*
    copy upper triangular matrix into in_buffer
    */
  
  d_ptr = in_buffer;
  for ( jndx = 0; jndx < nvecsL; jndx++ ){
    ii = mapvecL[jndx];
    ii = ll - ii;
    dcopy_ ( &ii, colL[jndx], &IONE, d_ptr, &IONE);
    d_ptr += ii;
  }
  
#ifdef DEBUG1
  fprintf (stderr, "8 after maxsz me = %d, i = %d, n=%d, m=%d \n", me, i, *n, *m);
#endif
  
  /*
    compute local matrix multiply
    */
  
  /*
    scratch stores the local copy of the full matrix
    */
  
  
  d_ptr = scratch;
  
#ifdef DEBUG
  for ( k = 0; k < nvecsF; k++ )
    for ( i = 0; i < *n; i++ )
      fprintf (stderr, "before mxm_ll1 me = %d F[%d][%d] = %f \n", me, k, i, colF[k][i]);
  
  
  for ( k = 0; k < nvecsF; k++ )
    for ( i = 0; i < *n; i++ )
      fprintf (stderr, "before mxm_ll1 me = %d F[%d][%d] = %f \n", me, k, i, *(d_ptr++));
  
  for ( k = 0; k < nvecsL; k++ )
    for ( i = 0; i < *n - mapvecL[k]; i++ )
      fprintf (stderr, "before mxm_ll1 me = %d L[%d][%d] = %f \n", me, k, i, colL[k][i]);
#endif
  
  F_ptr = (DoublePrecision *) scratch;
  
  for ( k = 0; k < nvecsF; k++ ){
    L_ptr = in_buffer;
    fff = mapvecF[k];
    for ( jndx = 0; jndx < nvecsL; jndx++ ) {
      i = mapvecL[jndx];
      if ( i >= fff ) {
	isize = ll - i;
	t = *(F_ptr + i);
	daxpy_ ( &isize, &t, L_ptr, &IONE, colF[k] + i , &IONE);
      }
      else {
	isize = ll - fff;
	t = *(F_ptr + i);
	daxpy_ ( &isize, &t, L_ptr+fff-i, &IONE, colF[k] + fff , &IONE);
      }
      L_ptr += (ll - i);
    }
    F_ptr += ll;
  }
  
#ifdef DEBUG1
  fprintf (stderr, "9 after maxsz me = %d, i = %d, n=%d, m=%d nproc = %d  \n", me, i, *n, *m, nproc);
  for (k = 0; k < nproc; k++ )
    fprintf(stderr, " me = %d proclist[%d] = %d \n", me, k, proclist[k]);
#endif
  
  
  /*
    buffer     = out_buffer + maxsz;
    */
  
  isize = ci_size_ ( &me, &ll, mapL);   /* the number Q vectors that I own */
  osize = isize;

  me_indx += nproc;      /* adjusting for mod */
  
  for ( i = 0; i < nproc-1 ; i++ ) {
    dcopy_ (&osize, in_buffer, &IONE, out_buffer, &IONE);
    if (( me_indx  % 2 ) ==  0 ) {
      
      rsize = osize * sizeof(DoublePrecision);
      if ( rsize != 0 ){ 
	rsize = mxwrit_ ( out_buffer, &rsize, &next_proc, &i );
      }
      
      indx = (( me_indx - i - 1) % nproc );

      indx = proclist[indx]; /* the current processor id number */
      
#ifdef DEBUG1
      fprintf (stderr, "91 after maxsz me = %d, i = %d, n=%d, m=%d nproc = %d indx = %d  \n", me, i, *n, *m, nproc, indx );
      for (k = 0; k < nproc; k++ )
	fprintf(stderr, " me = %d proclist[%d] = %d indx = %d \n", me, k, proclist[k], indx);
#endif

      nvecsL = fil_mapvec_ ( &indx, &ll, mapL, mapvec_in);
      isize = ci_size_ ( &indx, &ll, mapL);   /* the number Q vectors that I own */


      
      /*
	mapvec_in contains the mapvecQ for processor proclist[i])
	*/
      
      rsize = isize * sizeof(DoublePrecision);
      if ( rsize != 0 )
	rsize = mxread_( in_buffer, &rsize, &last_proc, &i );
    }
    else {
      indx = ( me_indx - i -1 ) % nproc;
      
      indx = proclist[indx]; /* the current processor id number */

#ifdef DEBUG1
      fprintf (stderr, "92 after maxsz me = %d, i = %d, n=%d, m=%d nproc = %d  \n", me, i, *n, *m, nproc);
      for (k = 0; k < nproc; k++ )
	fprintf(stderr, " me = %d proclist[%d] = %d indx = %d \n", me, k, proclist[k], indx);
#endif

      nvecsL = fil_mapvec_ ( &indx, &ll, mapL, mapvec_in);
      isize = ci_size_ ( &indx, &ll, mapL);   /* the number Q vectors that I own */


      
      /*
	mapvec_in contains the mapvecQ for processor proclist[i])
	*/
      
      rsize = isize * sizeof(DoublePrecision);
      
      if ( rsize != 0 )
	rsize = mxread_( in_buffer, &rsize, &last_proc, &i );
      
      rsize = osize * sizeof(DoublePrecision);
      if ( rsize != 0 )
	rsize = mxwrit_ ( out_buffer, &rsize, &next_proc, &i );
    }
    
    osize = isize;
    
    F_ptr = (DoublePrecision *) scratch;
    
    for ( k = 0; k < nvecsF; k++ ){
      L_ptr = in_buffer;
      fff = mapvecF[k];
      for ( jndx = 0; jndx < nvecsL; jndx++ ) {
	ii = mapvec_in[jndx];
	if ( ii >= fff ) {
	  isize = ll - ii;
	  t = *(F_ptr + ii);
	  daxpy_ ( &isize, &t, L_ptr, &IONE, colF[k] + ii , &IONE);
	}
	else {
	  isize = ll - fff;
	  t = *(F_ptr + ii);
	  daxpy_ ( &isize, &t, L_ptr+fff-ii, &IONE, colF[k] + fff , &IONE);
	}
	L_ptr += (ll - ii);
      }
      F_ptr += ll;
    }
  }
  
#ifdef DEBUG1
  fprintf(stderr, " Out of mxm_llx *********************** %d \n", me );
#endif
  
  return;
}
