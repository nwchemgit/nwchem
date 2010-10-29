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

void mxm5x( n, rowU, mapU, m, colF, mapF, iscratch, scratch)
     Integer *n, *mapU, *mapF,  *m, *iscratch;
     DoublePrecision **colF, **rowU, *scratch;
{
  
/**************************************************************
 *
 * Subroutine mxm5x
 *
   This subroutine does the matrix multiplication F <- U * F
       
   where U is a n x n upper trianguler matrix in packed storage format,
   and row distributed,
   and matrix F is a general n x m matrix distributed by columns.
       
   ARGUMENTS
   ---------
   In the following:

     INTEGER          = "pointer to Integer"
     DOUBLE PRECISION = "pointer to DoublePrecision"

     me     = this processor's id (= mxmynd_())
     nprocs = number of allocated processors ( = mxnprc_())
     nvecsU = number of entries in mapU equal to me
              (= count_list( me, mapU, n ))
     neleU_max = maximum number of elements of U owned by any
                 processor.
                 (= max over i of ci_size_( i, n, mapU) )
     nvecsF = number of entries in mapF equal to me
              (= count_list( me, mapF, m ))
     sDP    = sizeof( DoublePrecision )

   n ...... (input) INTEGER
            The number of rows and columns of U, and the number
            of rows of F.

   rowU ... (input) array of pointers to DoublePrecision,
                    length (nvecsU)
            The part of matrix U owned by this processer stored
            in packed format, i.e., rowU[i] points to the diagonal
            element of the i-th row of U
            owned by this processor, i = 0 to nvecsU-1.
                
   mapU ... (input) INTEGER array, length (n)
            The i-th row of U is 
            owned by processor mapU[i], i = 0 to n-1.

   m ...... (input) INTEGER
            m is the number of columns in F.
       
   colF ... (input/output) array of pointers to DoublePrecision,
                           length (nvecsF)

            On Entry:
              The part of matrix F owned by this processer stored
              in packed format by columns, i.e., colF[i] points
              to the first element of the i-th column of F
              owned by this processor, i = 0 to nvecsF-1.
                
            On Exit:
              The result matrix U * F stored in the same manner as
              F on entry.

   mapF ... (input) INTEGER array, length (m)
            The i-th column of F is 
            owned by processor mapF[i], i = 0 to m-1.
                
   iscratch .. (workspace) INTEGER array, length ( 3*n+2*m +1 )

   scratch.... (workspace) DOUBLE PRECISION array, 
                        length ( maximum ( n+m, mxlbuf_()/sDP + 1,
                                           (nvecsF * n + 2 * neleU_max )
*/
  static Integer ONE = 1;
  Integer jndx;
  Integer i, k, me, isize, *iptr, indx;
  Integer nvecsU, nvecsF, ii;
  Integer j, linfo, ll, mm, nproc;
  Integer osize, rsize, idummy;		/* */
  
  Integer *mapvecF, *mapvecU;
  Integer *iscrat, *proclist;
  Integer *mapvec_in, me_indx, last_proc, next_proc, maxsz;
  
  DoublePrecision *buffer, *d_ptr, *in_buffer, *out_buffer;
  DoublePrecision *F_ptr, *U_ptr, t;

  
  
  /*
    blas calls
    */
  
  extern void dcopy_(), daxpy_();
  extern DoublePrecision ddot_();

  extern Integer mxmynd_();
  extern void xerbla_();
  extern void mapdif1_();
  extern void g_exit_();
  extern void pxerbla2_();
  extern Integer fil_mapvec_();
  extern void chol_pipe_bcast();
  extern void gshellsort_();
  extern Integer indxL();
  extern void fil_dbl_lst();

  
  /* 
    mxsubs calls
    */
  
  extern Integer mxwrit_(),  mxread_();
  extern Integer count_list();
  extern Integer fil_mapvec_();
  extern Integer ci_size_();
  extern Integer reduce_list2();
  extern Integer peigs_cmod_();
  
  i = 0;
  me = mxmynd_();

#ifdef DEBUG1
  fprintf(stderr, " Entering mxm5x ********************* me = %d \n ", me );
#endif
  
  if ( m == NULL )
    i = -4;
  else
    if ( n == NULL ) {
      i = -1;
      xerbla_( "MXM5 ", &i);
    }
    else if ( *n < 1 )
      i = -1;
    else {
      iscrat = mapU;
      for ( j = 0; j < *n; j++ ) {
	if ( iscrat == NULL ) {
	  i = -3;
	  xerbla_( "MXM5 \n", &i);
	}
	else
	  iscrat++;
      }
      iscrat = mapF;
      for ( j = 0; j < *m; j++ ) {
	if ( iscrat == NULL ) {
	  i = -6;
	  xerbla_( "MXM5 \n", &i);
	}
	else
	  iscrat++;
      }
    }
  
  /*
    at this point inputs are minimally acceptable
    
    check to see if mapU and mapF are the same set of processors
    */
  
  iscrat = iscratch;
  mapdif1_( n, mapU, m, mapF, iscrat, &j );
  
  /*
   * if ( j != 0 ) {
   *   i = -3;
   * }
   */
  
#ifdef DEBUG1
  fprintf (stderr, "me = %d, i = %d, n=%d, m=%d \n", me, i, *n, *m);
  for ( j = 0; j < *n; j++ )
    fprintf (stderr, "me = %d, mapU[%d]= %d \n", me, j, mapU[j] );
  for ( j = 0; j < *m; j++ )
    fprintf (stderr, "me = %d, mapF[%d]=%d, \n", me, j, mapF[j] );
#endif
  
  me = mxmynd_();
  nvecsU = count_list ( me, mapU, n);
  nvecsF = count_list ( me, mapF, m);
  
  for ( j = 0; j < nvecsU; j++ )
    if ( rowU[j] == NULL ) {
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

    g_exit_( &i, "Mapping problem or memory assignment problem MXM5X \n", mapU, n, iscratch, scratch );
    g_exit_( &i, "Mapping problem or memory assignment problem MXM5 \n", mapF, n, iscratch, scratch );
    
    
    linfo = 0;
    ll = *n * sizeof(Integer);
    pxerbla2_( &ll, (char *) mapU, mapU, n, iscrat, &i );
    linfo = min(linfo, i);
    ll = *m * sizeof(Integer);
    
    pxerbla2_( &ll, (char *) mapF, mapU, n, iscrat, &i);
    
    linfo = min(linfo, i);
    
    g_exit_( &linfo, "Mapping inconsistancies calling MXM5X \n", mapU, n, iscratch, scratch );
    g_exit_( &linfo, "Mapping inconsistancies calling MXM5X \n", mapF, n, iscratch, scratch );
    */
  
  me = mxmynd_();
  
  ll = *n;
  mm = *m; 
  
  iscrat = iscratch;
  mapvecU = iscrat;

#ifdef DEBUG
  fprintf(stderr, " in mxm5x me = %d \n", me);
#endif
  
  nvecsU = fil_mapvec_( &me, &ll, mapU, mapvecU );
  iscrat += nvecsU;
  mapvecF = iscrat;
  nvecsF = fil_mapvec_( &me, m, mapF, mapvecF );
  iscrat += nvecsF;
  
#ifdef DEBUG1
  k = -1;
  for ( j = 0; j < *n; j++ ){
    if( mapU[j] == me ) {
      k++;
      for ( i = 0; i < *n-j; i++ )
        fprintf (stderr, "me = %d, rowU[%d][%d] = %g \n", me, k, i, rowU[k][i] );
    }
  }
  for ( j = 0; j < nvecsF; j++ )
    for ( i = 0; i < *n; i++ )
      fprintf (stderr, "me = %d, colF[%d][%d] = %g \n", me, j, i, colF[j][i] );
#endif


  if ( nvecsU + nvecsF == 0 )
    return;
  
  iptr = (Integer *) scratch;
  for ( i = 0; i < ll; i++ )
    *(iptr++) = mapU[i];
  for ( i = 0; i < mm; i++ )
    *(iptr++) = mapF[i];
  
  iptr = (Integer *) scratch;
  
  /*
    nproc = reduce_list( ll+mm, iptr, iscrat);
    */
  
  nproc = reduce_list2( ll+mm, iptr, iscrat);
  proclist = iscrat;
  iscrat += nproc;
  mapvec_in = iscrat;
  
  gshellsort_( &nproc, proclist);
  
  indx = indxL ( me, nproc, proclist);
  me_indx = indx;
  indx += nproc;
  
  /*
    last_proc = (indx - 1) % nproc;
    */
  
  idummy = indx - 1;
  idummy += nproc;
  last_proc = (idummy+nproc) % nproc;
  last_proc = proclist[last_proc];
  
  /*
    next_proc = (indx + 1) % nproc;
    */
  
  idummy = indx + 1;
  idummy += nproc;
  next_proc = (idummy + nproc ) % nproc;
  next_proc = proclist[next_proc];
  
  maxsz = 0;
  isize = 0;
  for ( indx = 0; indx < nproc; indx++ ) {
    isize = ci_size_( &proclist[indx], &ll, mapU);   /* the number Q vectors that I own */
    maxsz = max( isize, maxsz );
  }
  
  buffer = (DoublePrecision *) scratch;
  indx = 0;
  for ( jndx = 0; jndx < nvecsF; jndx++ ){
    dcopy_( &ll, colF[jndx], &ONE, &buffer[indx], &ONE);
    indx += ll;
  }
  buffer += indx;
  
  for ( jndx = 0; jndx < nvecsF; jndx++ ){
    d_ptr = colF[jndx];
    fil_dbl_lst ( ll, d_ptr, 0.0e0);
  }
  
  in_buffer = buffer;
  out_buffer = buffer + maxsz;
  
  /*
    copy upper triangular matrix into memory
    */
  
  d_ptr = in_buffer;
  for ( jndx = 0; jndx < nvecsU; jndx++ ){
    ii = mapvecU[jndx];
    ii = ll - ii;
    dcopy_( &ii, rowU[jndx], &ONE, d_ptr, &ONE);
    d_ptr += ii;
  }
  
  /*
    compute local matrix multiply
    */
  
  F_ptr = scratch;
  for ( k = 0; k < nvecsF; k++ ){  
    U_ptr = in_buffer;
    for ( jndx = 0; jndx < nvecsU; jndx++ ) {
      i = mapvecU[jndx];
      isize = ll - i;
      t = ddot_( &isize, &F_ptr[i], &ONE, U_ptr , &ONE);
      colF[k][i] = t;
      U_ptr += isize;
    }
    F_ptr += ll;
  }
  
  buffer     = out_buffer + maxsz;
  
  isize = ci_size_( &me, &ll, mapU);   /* the number Q vectors that I own */
  osize = isize;
  
  me_indx += nproc;
  
  for ( i = 0; i < nproc-1 ; i++ ) {
    dcopy_(&osize, in_buffer, &ONE, out_buffer, &ONE);
    idummy = 2;
    if ( (me_indx + nproc) % idummy ==  0 ) {
      rsize = osize * sizeof(DoublePrecision);
      if ( rsize != 0 )
	rsize = mxwrit_( out_buffer, &rsize, &next_proc, &i );
      
      /*
        indx = ( me_indx -i -1) % nproc;
        */
      idummy = me_indx - i - 1;
      idummy += nproc;
      indx = (idummy + nproc ) % nproc;
      indx = proclist[indx]; /* the current processor id number */
      nvecsU = fil_mapvec_( &indx, &ll, mapU, mapvec_in);
      isize = ci_size_( &indx, &ll, mapU);   /* the number Q vectors that I own */
      
      /*
	mapvec_in contains the mapvecQ for processor proclist[i])
	*/
      
      rsize = isize * sizeof(DoublePrecision);
      if ( rsize != 0 )
	rsize = mxread_( in_buffer, &rsize, &last_proc, &i );
    }
    else {
      /*
        indx = ( me_indx - i -1 ) % nproc;
        */
      idummy = me_indx - i - 1;
      idummy += nproc;
      indx = (idummy + nproc) % nproc;
      indx = proclist[indx]; /* the current processor id number */
      nvecsU = fil_mapvec_( &indx, &ll, mapU, mapvec_in);
      isize = ci_size_( &indx, &ll, mapU);   /* the number Q vectors that I own */
      
      /*
	mapvec_in contains the mapvecQ for processor proclist[i])
	*/
      
      rsize = isize * sizeof(DoublePrecision);

      if ( rsize != 0 )
	rsize = mxread_( in_buffer, &rsize, &last_proc, &i );
      
      rsize = osize * sizeof(DoublePrecision);
      if ( rsize != 0 )
	rsize = mxwrit_( out_buffer, &rsize, &next_proc, &i );
    }
    
    osize = isize;
    
    F_ptr = scratch;
    for ( k = 0; k < nvecsF; k++ ){  
      U_ptr = in_buffer;
      for ( jndx = 0; jndx < nvecsU; jndx++ ) {
	ii = mapvec_in[jndx];
	isize = ll - ii;
	t = ddot_( &isize, &F_ptr[ii], &ONE, U_ptr , &ONE);
	colF[k][ii] = t;
	U_ptr += isize;
      }
      F_ptr += ll;
    }
  }

  return;
}

