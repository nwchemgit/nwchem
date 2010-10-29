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
/* **************************************************
 *
 * PeIGS routine inverseL
 * the auxiliary communication routine 
 *        
 *        void pipe_bcst_prev_col(buf, len, type, root, c_indx, map, scratch)
 * and 
 *        Integer prev_column( node, map, indx)
 *
 * also live here
 *
 *  this routine computes the inverse of a real lower triangular
 *  matrix stored in the packed storage format
 *
 *  on exit, it overwrites the current lower triangular matrix with the the inverse
 * of the lower triangular matrix
 *
 *
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "globalp.c.h"

#define min(a,b) ((a) < (b) ? (a) : (b))

#define C_sign(a, b) ((b) < 0.0 ? -(a) : (a))
#define i_have(i, numprocs, ime) ( (i) % numprocs == ime)

void pipe_bcst_prev_col(buf, len, type, root, c_indx, map, scratch)
     char *buf;
     Integer len;
     Integer type;
     Integer root, c_indx;
     Integer *map, *scratch;
     /* 
       In this version of broadcast, we are broadcasting to a list of
       processors holding columns with column numbers less than c_indx
       
       Duplications in the processor list is removed
       
       c_indx  = current index; the array map is assumed to have size at least
       
       map     = is a ptr to the beginning of a
                 list of locations of the columns/rows
                 of the matrices
       
       root    = the node that started it all
       
       scratch = scratch integer memory of length at least map
       
       
       */
{
  Integer num_indx, me;
  Integer *junk, *ptr, k, isize;
  extern Integer mxmync_();
  extern Integer mxmynd_();
  extern Integer mxnprc_();
  extern Integer mxwrit_();
  extern Integer mxread_();
  

  extern Integer in_list();
  
  Integer prev_column();
  
  /*
    reduce the list of processors to a minimum list holding previous column numbers
    */
  
  me = mxmynd_();
  
  if (  mxnprc_() == 1 )
    return;
  
  if ( c_indx == 0 )
    return;
  
  junk = scratch;
  *junk = root;
  num_indx = 1;
  
  /*
     reduce the list to only
     */
  
  for ( k = 0; k < c_indx; k++ ) {
    if (( map[k] != root ) && ( in_list( &map[k], junk, &num_indx) == 0 )) {
      junk[num_indx] = map[k];
      num_indx++;
    }
  }
  
  /*
     num_indx holds the number of processors holding previous columns
     */
  
  /*
     junk[1:num_indx] now contains the list of processors in order without root and duplication
     in sequence of the method they appear according to map[*]
     */
  
  
  if ( num_indx == 1 )
    return;  /* root is the only node */
  
  ptr = junk + 1;
  if ( me == root )
    isize = mxwrit_( buf, &len, ptr, &type);
  
  for ( k = 1; k < num_indx; k++ ){
    if ( me == junk[k] )	{
      isize = mxread_( buf, &len, &junk[k-1] , &type);
      if ( k < ( num_indx - 1 ) ) {
	isize = mxwrit_( buf, &len, &junk[k+1] , &type);
      }
    }
  }
  return;
}

/*
  
  This routine inverts a lower triangular matrix using a list-column distributed
  storage format.  On output the inverse matrix overwrites the input matrix
  with the same storage mapping.
  
  */

void inverseL ( msize, col, map, iwork, buffer, info)
     Integer *msize, *map, *iwork, *info;
     DoublePrecision **col, *buffer;
     
     /*

       This routine computes the inverse of a lower triangular matrix.
       On output, the lower triangular matrix is overwritten by its inverse.
       The matrix is assume to be column wrapped and distributed according
       to the array map.  Thus, map[i] = processor id which holds column
       i of the matrix.

       Arguments:
       
       msize = (input) integer, size of the matrix
       
       col = (input/1-D array of pointers to DoublePrecision precision numbers)
             DoublePrecision pointer to a DoublePrecision array of matrices
	     col[j] -> the address of the j-th vector that this processor owns
	     
       map = (input/integer array) integer array, length n, such that for the i-th column,
	     map[i] = real processor id holding this column
       
       iscratch = (scratch/integer) integer scratch space of size 3 * msize + 3

       buffer = (scratch/DoublePrecision precision) DoublePrecision precision array of
                length MAX{ msize + 1, mxlbuf_() / sizeof(DoublePrecision) + 1 }
       
       info  = (output/integer) if info == 0 normal
                                if info == i > 0 then the i-th diagonal element
				( in Fortran notation ) is smaller than the dlamch("s"), the
				safe inverse limit: overflow may be occurred.
				
     *
     */
{
  static Integer IONE = 1;
  
  Integer n, i, j, me, isize, iL;
  Integer linfo, nvecs;
  Integer *myvecs, *iscrat;
  
  DoublePrecision temp, *q, safeulp;
  extern Integer mxmynd_();
  extern void pipe_bcst_prev_col();
  extern void xerbla_();
  extern void bbcast00();
  extern Integer count_list();
  extern void pxerbla2_();
  extern void g_exit_();
  extern Integer fil_mapvec_();
  extern void dcopy_();
  extern void dscal_();
  extern void daxpy_();


  Integer prev_column();
  Integer num_procs;
  extern Integer reduce_list(), reduce_list2();
  Integer length;

/*
  extern void iCC_work_init(), iCC_bcast();
*/
  
  i = 0;
  me = mxmynd_();
  
  if ( msize == NULL ) {
    i = -1;
    xerbla_( "INVERSE ", &i);
  }
  else
    if ( *msize < 1 ){
      i = -1;
      xerbla_( "INVERSE ", &i);
    }
    else {
      iscrat = map;
      for ( j = 0; j < *msize; j++ ) {
	if ( iscrat == NULL ) {
	  i = -3;
	  xerbla_( "INVERSE \n", &i);
	}
	else
	  iscrat++;
      }
    }
  
  
  /*
    at this point inputs are minimally acceptable
    */
  
  nvecs = count_list ( me, map, msize);
  
  for ( j = 0; j < nvecs; j++ ) {
    if ( col[j] == NULL ) {
      linfo = -3;
      i = min(linfo, i);
      break;
    }
  }
  
  iscrat = iwork;
  *iscrat = *msize;
  
  n = *msize;  ;     /* set dim of the matrix */
  for ( j = 0; j < n; j++ )
    iscrat[j+1] = map[j];
  
  linfo = 0;
  j = n + 1;
  j = j*sizeof(Integer);
  
  pxerbla2_( &j, (char *) iscrat, map, msize, iwork + n + 1, &linfo );
  
  linfo = -abs(linfo);
  
  g_exit_( &linfo, "Mapping problem or memory assignment problem INVERSE \n", map, &n, iwork, buffer );
  
  me = mxmynd_();
  iscrat = iwork;
  myvecs = iscrat;
  nvecs = fil_mapvec_( &me, &n, map, myvecs );
  iscrat += nvecs;
  
  if ( nvecs == 0 )
    return;
  
  
  safeulp = DLAMCHS;
  
  *info = 0;
  
  for ( i = 0; i < nvecs; i++ ) {
    if ( *col[i] < safeulp )
      *info = i+1;
    col[i][0] = 1.0e0/col[i][0];
  }
  
  num_procs = reduce_list2( n, map, iscrat );
  
#ifdef MESH
  i_random(num_procs, iscrat, iscrat+num_procs+1);
#endif
  
  /*
    used now for broadcast
    */
  
  
  iL = 0;
  length = 2*n;
  
  /*
    iCC_work_init((unsigned long) length );
    */
  
  for ( i = 0; i < n; i++ ){
    if ( map[i] == me ){
      isize = n-i;
      dcopy_( &isize, col[iL], &IONE, buffer, &IONE);
      bbcast00( (char *) buffer, (n-i)*sizeof(DoublePrecision), i*100, me, num_procs, iscrat );
      temp = -col[iL][0];
      isize = n-i-1;
      dscal_( &isize, &temp, &col[iL][1], &IONE);
      iL++;
    }
    
    if ( map[i] != me) 
      bbcast00( (char *) buffer, (n-i)*sizeof(DoublePrecision), i*100, map[i], num_procs, iscrat );

    
    /*
      the buffer holds ( 1/l[i,i], l[i+1, i], l[i+2, i], ..., l[n-1, i])
      */
    
    isize = n-i-1;      
    for ( j = 0; j < nvecs; j++ ){
      if ( myvecs[j] < i ){
	q = col[j] + i - myvecs[j];
	*q  *= *buffer;
	temp = - *q;
	q++;
	daxpy_( &isize, &temp, &buffer[1], &IONE, q, &IONE);
      }
    }
  }
  return;
}

