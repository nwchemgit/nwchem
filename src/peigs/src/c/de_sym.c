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
#define max(a,b) ((a) > (b) ? (a) : (b))
#define min(a,b) ((a) < (b) ? (a) : (b))

/*
  Internal Peigs subroutine:
  
  void de_sym(n, msgtype, buffer, matrix, map, matrix2, iwork, buf_ptr )
  Integer n, msgtype, *map, *iwork;
  DoublePrecision *buffer, **matrix, **matrix2, **buf_ptr;
  
  Purpose: to duplicate the lower half of a symmetric matrix to 
  the upper half of a symmetric matrix.  The matrix is assumed to
  be stored using a column distribution across a list of processors given
  in map.
  
  */

void de_sym(n, msgtype, buffer, matrix, map, matrix2, iwork, buf_ptr )
     Integer n, msgtype, map[], iwork[];
     DoublePrecision buffer[], **matrix, **matrix2, **buf_ptr;
{
  /*
    
    This routine duplicates the lower half of a symmetric matrix
    to have the upper half:de-symmetrize the packed storage format.
    This part is similiar to half of a transpose.
    
    dimension buf(lenbuf, 0:*), buffer(lenbuf)
    integer lengthbuf
    
    Every process needs to exchange buffers with every other process
    
    buf(lenbuf, 0:nproc-1) - on input contains buffer to be sent to
    process iproc=0:nproc-1
    - on output contains buffer received from
    process iproc=0:nproc-1
    lenbuf - length of each buffer in REAL*8 words
    buffer - workspace of length lenbuf
    
    Assume processes are directly connected ... sequence of
    exchanges is designed for maximum overlap of messages so
    will complete in time O(nproc*lenbuf) if the hardware permits
    (and all is done asynand messages long enuf to overcome
    startup overhead ... get real!).
    
    This may need modification to reflect specifihardware toplogy.
    e.g. for 6 processors have 5 steps each with 3 exchanges in parallel
    
    0-1     0-2     0-3     0-4     0-5
    5-2     1-3     2-4     3-5     4-1
    4-3     5-4     1-5     2-1     3-2
    
    */
  
  Integer nvecs, *iscrat, ip, ig, lenb, i, me_indx, idest;
  Integer me, nprocs, ngroup, nvecs_out;
  Integer *mapvec, *map_out, *proclist;
  Integer load_up(), un_load(), un_load_size();
  DoublePrecision *read_buffer, *write_buffer;
  
  extern Integer mxnprc_ (), mxmynd_ ();
  extern void mxread_ (), mxwrit_();
  
  extern void pairup_ ();
  extern Integer fil_mapvec_ ();
  extern Integer reduce_list2();
  extern Integer indxL();
  extern void pairup_ ();
  extern Integer ci_size_ ();

  me = mxmynd_();

  iscrat = iwork;
  proclist = iscrat;
  nprocs = reduce_list2( n, map, iscrat);
  iscrat += nprocs;
  me_indx = indxL( me, nprocs, proclist);
  
  ngroup =(nprocs % 2);
  ngroup += nprocs;
  
#ifdef DEBUG1
  fprintf(stdout, " de_sym me = %d me_indx %d ngroup %d n %d nprocs %d   \n", me, me_indx, ngroup, n, nprocs);
  fflush(stdout);
#endif
  
  /*
   */
  
  mapvec = iscrat;
  nvecs = fil_mapvec_ ( &me, &n, map, mapvec );
  iscrat += nvecs;
  
  if ( nvecs == 0 )
    return;
  
#ifdef DEBUG2  
  fprintf(stderr, " de_sym 2 %d \n", me);
  
  for ( ig = 0; ig < nvecs; ig++ )
    for ( ip = 0; ip < n - mapvec[ig]; ip++ ) {
      fprintf(stderr, " me = %d matrix [ %d ] [ %d ] = %g \n", me, ig, ip, matrix[ig][ip]);
    }
  
  for ( ig = 0; ig < nvecs; ig++ ) {
    for ( ip = 0; ip < n; ip++ ) {
      fprintf(stderr, " me = %d matrix2 [ %d ] [ %d ] = %g \n", me, ig, ip, matrix2[ig][ip]);
    }
  }
#endif
  
  i = 0;
  ig = 0;
  for ( ip = 0; ip < nprocs; ip++ ) {
    ig = ci_size_( &iwork[ip], &n, map);
    i = max( ig, i );
  }
  
  read_buffer  =  buffer;
  write_buffer =  &read_buffer[i];
  
  lenb = load_up(nvecs, mapvec, nvecs, mapvec, matrix, write_buffer);
  
  un_load(nvecs, mapvec, nvecs, mapvec, matrix2, write_buffer);
  
  ngroup = nprocs;
  ngroup += (nprocs % 2);

  map_out = iscrat;
  for ( ig = 1; ig <=  (ngroup-1); ig++ ) {
    pairup_(&ngroup, &me_indx, &ig, &ip);
    fflush(stdout);
    
    if (ip <  nprocs) {
      if (ip > me_indx) {
	idest = proclist[ip];
	nvecs_out = fil_mapvec_( &idest, &n, map, map_out );
	lenb = load_up(nvecs, mapvec, nvecs_out, map_out, matrix, write_buffer);
	lenb *= sizeof(DoublePrecision);
	if ( lenb != 0 ) {
	  mxwrit_ (write_buffer, &lenb, &idest, &msgtype);
	}
	
	/*
	 */
	
	lenb =un_load_size(nvecs, mapvec, nvecs_out, map_out);
	lenb *= sizeof(DoublePrecision);	
	if ( lenb != 0 ) {
	  mxread_ (read_buffer, &lenb, &idest, &msgtype);
	  un_load(nvecs, mapvec, nvecs_out, map_out, matrix2, read_buffer);
	}
      }
      else {
	idest = proclist[ip];
	nvecs_out = fil_mapvec_ ( &idest, &n, map, map_out );
	lenb =un_load_size(nvecs, mapvec, nvecs_out, map_out);
	lenb *= sizeof(DoublePrecision);	
	if ( lenb != 0 )  {
	  mxread_ (read_buffer, &lenb, &idest, &msgtype);
	  un_load(nvecs, mapvec, nvecs_out, map_out, matrix2, read_buffer);
	}
	
	lenb = load_up(nvecs, mapvec, nvecs_out, map_out, matrix, write_buffer);
	lenb *= sizeof(DoublePrecision);
	if ( lenb != 0 ){
	  mxwrit_ (write_buffer, &lenb, &idest, &msgtype);
	}
      }
    }
  }
  
  
#ifdef DEBUG  
  for ( ig = 0; ig < nvecs; ig++ ) {
    for ( ip = 0; ip < n; ip++ ) {
      fprintf(stderr, " me = %d matrix [ %d ] [ %d ] = %g \n", me, ig, ip, matrix2[ig][ip]);
    }
  }
  
  fprintf(stderr, " ******************************* exiting de_sym me = %d \n", me );
#endif
  return;
}

Integer load_up(nvecs, mapvec, nvecs_out, map_out, matrix, buffer)
     Integer nvecs, *mapvec, nvecs_out, *map_out;
     DoublePrecision **matrix, *buffer;
{
  Integer indx, jndx, j, k, icount;
  DoublePrecision *dptr, *column_ptr;
  
  dptr = buffer;
  icount = 0;
  for ( indx = 0; indx < nvecs; indx++ ) {
    k = mapvec[indx];
    column_ptr = matrix[indx];
    for ( jndx = 0; jndx < nvecs_out; jndx++ ) {
      j = map_out[jndx];
      if ( j > k ){
	*dptr = column_ptr[j-k];
	dptr++;
	icount++;
      }
    }
  }
  return(icount);
}

Integer un_load(nvecs, mapvec, nvecs_in, map_in, matrix, buffer)
     Integer nvecs, *mapvec, nvecs_in, *map_in;
     DoublePrecision **matrix, *buffer;
{
  Integer indx, jndx, j, k, icount;
  DoublePrecision *dptr;
  
  dptr = buffer;
  icount = 0;
  
  for ( jndx = 0; jndx < nvecs_in; jndx++ ) {
    j = map_in[jndx];
    for ( indx = 0; indx < nvecs; indx++ ) {
      k = mapvec[indx];
      if ( j < k ){
	/*
	  fprintf(stderr, " storing me = %d matrix[ %d] [%d]  = %f \n", mxmynd_ (), indx, j , *dptr);
	  */
	matrix[indx][j] = *dptr;
	dptr++;
      }
    }
  }
  return(0);
}

Integer un_load_size(nvecs, mapvec, nvecs_in, map_in)
     Integer nvecs, *mapvec, nvecs_in, *map_in;
{
  Integer indx, jndx, j, k, icount = 0;
  
  for ( indx = 0; indx < nvecs; indx++ ) {
    k = mapvec[indx];
    for ( jndx = 0; jndx < nvecs_in; jndx++ ) {
      j = map_in[jndx];
      if ( j < k )
	icount++;
    }
  }
  return(icount);
}

