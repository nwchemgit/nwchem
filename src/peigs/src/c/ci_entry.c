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

/* **************************************************************************
 *
 *	 PeIGS Utility Routines
 *
 * fortran and C callable
 *
 this file contains subroutines:

 ci_entry :returns index a[j][i] for a 1-D data array
 storing the lower triangular part of a
 symmetric matrix
 
 ci_entry_ :returns index a(i,j) in Fortran for a 1-D data array
            storing the lower triangular part of a
            symmetric matrix
 
 ci_size_ :  returns the total memory required
             for "standard" column
              wrapping of a symmetric matrix
 
fil_mapvec_ : from map construct a correspondance
              mapvec: the real column/row index that the matrix
              on this processor owns

*/


Integer ci_entry (me, n, i, j, map)
     Integer *me, *n, *i, *j, *map;
{
  Integer indx, jndx, label;
  Integer *iptr;
  extern void xerbla_();
  
  /*
    PeIGS utility routine
    
    this routine returns the C index of the a(i,j)
    for a symmetric matrix with i >= j
    on processor me stored using a 1-D array
    
    if (processor = me)
    doesn't own this element -1 is returned
    
    Argument:
    
    me        = (input/integer) processor id
    
    n          = (input/integer) dimension of the matrix
    
    i           = (input/integer) the row index of the element a[j][i]
    
    j           = (input/integer) the column index of the element a[j][i]
    
    map      = (input/integer array ) integer array of length n
    map[j] = the processor id holding column j
    
    */
  
  if ( me == NULL ) {
    indx = -1;
    xerbla_("ci_entry\n", &indx);
  }
  
  if ( n == NULL ) {
    indx = -2;
    xerbla_("ci_entry\n", &indx);
  }
  else if ( *n < 1 ) {
    indx = -2;
    xerbla_("ci_entry\n", &indx);
  }
  
  if ( i == NULL ) {
    indx = -3;
    xerbla_("ci_entry\n", &indx);
  }
  
  if ( j == NULL ) {
    indx = -4;
    xerbla_("ci_entry\n", &indx);
  }
  
  if (map == NULL) {
    printf(" Peigs: ci_entry_  node %d : Mapping problem\n", *me);
    indx = -5;
    xerbla_("ci_entry\n", &indx);
  }
  
  jndx = -5;
  iptr = map;
  for ( indx = 0; indx < *j; indx++ ) {
    if ( (iptr++) == NULL )
      xerbla_("ci_entry\n", &jndx);
  }
  
  if ( *i < 0 ) {
    indx = -4;
    xerbla_("ci_entry\n", &indx);
  }
  
  if ( *j  < 0 ) {
    indx = -5;
    xerbla_("ci_entry\n", &indx);
  }
  
  if ( *i < *j ) {
    printf("PeIGS: ci_entry_ :node %d : i < j \n", *me);
    indx = -4;
    xerbla_("ci_entry\n", &indx);
  }
  
  jndx = 0;
  label = -1;
  
  for ( indx = 0; indx < *j; indx++ ) {
    if ( map[indx] == *me ) {
      jndx += *n - indx;
      label = 1;
    }
  }
  
  indx = *j;
  if ( map[ indx ] == *me ) {
    jndx += *i - indx;
    label = 1;
  }
  
  if ( label == -1 )
    return(-1);
  else
    return(jndx);
  
}

Integer ci_entry_(me, n, i, j, map)
     Integer *me, *n, *i, *j, *map;
{
  Integer indx, jndx, label, idummy, idummy2;
  Integer *iptr;
  extern void xerbla_();
  extern Integer mxmynd_();
  
  /*
    PeIGS utility routine
    
    this routine returns the Fortran index of the a(i,j)
    for a symmetric matrix with i >= j
    on processor me
    
    if (processor = me)
    doesn't own this element -1 is returned
    
    Argument:
    
    me        = (input/integer) processor id
    
    n          = (input/integer) dimension of the matrix
    
    i           = (input/integer) the row index of the element a[j][i]
    
    j           = (input/integer) the column index of the element a[j][i]
    
    map      = (input/integer array ) integer array of length n
    map[j] = the processor id holding column j
    
    */
  
  if ( me == NULL ) {
    indx = -1;
    xerbla_("ci_entry\n", &indx);
  }
  
  if ( n == NULL ) {
    indx = -2;
    xerbla_("ci_entry\n", &indx);
  }
  else if ( *n < 1 ) {
    indx = -2;
    xerbla_("ci_entry\n", &indx);
  }
  
  if ( i == NULL ) {
    indx = -3;
    xerbla_("ci_entry\n", &indx);
  }
  
  if ( j == NULL ) {
    indx = -4;
    xerbla_("ci_entry\n", &indx);
  }

  if (map == NULL) {
    printf(" Peigs: ci_entry_  node %d : Mapping problem\n", *me);
    indx = -5;
    xerbla_("ci_entry\n", &indx);
  }
  
  jndx = -5;
  iptr = map;
  for ( indx = 0; indx < *j; indx++ ) {
    if ( (iptr++) == NULL )
      xerbla_("ci_entry\n", &jndx);
  }
  
  if ( *i < 0 ) {
    indx = -4;
    xerbla_("ci_entry\n", &indx);
  }
  
  if ( *j  < 0 ) {
    indx = -5;
    xerbla_("ci_entry\n", &indx);
  }
  
  if ( *i < *j ) {
    printf( "PeIGS: ci_entry_ :node %d : i < j \n", *me);
    indx = -4;
    xerbla_("ci_entry\n", &indx);
  }
  
  jndx = 0;
  label = -1;
  iptr = map;
  idummy = *me;
  idummy2 = *n;
  
  for ( indx = 0; indx < *j; indx++ ) {
    if ( *(iptr++) == idummy ) {
      jndx += idummy2 - indx;
      label = 1;
    }
  }
  
  indx = *j;
  if ( map[ indx ] == *me ) {
    jndx += *i - indx;
    label = 1;
  }
  
  if ( label == -1 )
    return(-1);
  else
    return(jndx+1);
  
}

Integer ci_size_(me, n, map)
     Integer *me, *n, *map;
{
  /*
    Fortran and C callable
    
    returns the total number of DoublePrecision precision
    storage location required for storing the lower
    triangular part of a symmetric matrix distributed
    according the map
    
    me         = (input/integer) processor id
    
    n           = (input/integer) dimension of the matrix and also
    the length of the array map
    
    map(*)   = (input/integer array) length n array of processor ids
    
    */
  
  Integer indx, jndx, kndx;
  Integer *ptr, iam;
  
  extern void xerbla_();
  extern Integer mxnprc_();
  extern Integer mxmynd_();
  
  if ( me == NULL ) {
    indx = -1;
    xerbla_("ci_size\n", &indx);
  }

  iam = mxmynd_();
  
  if ( *me > mxnprc_() || *me < 0 ) {
    indx = -1;
    printf(" Node %d Error in ci_size arg 1 with me = %d mxnprc = %d \n", iam, *me, mxnprc_());
    xerbla_("ci_siz\n", &indx);
  }
  
  if ( n == NULL ) {
    indx = -2;
    xerbla_("ci_siz\n", &indx);
  }
  else if ( *n < 1 ) {
    indx = -2;
    xerbla_("ci_siz\n", &indx);
  }
  
  if (map == NULL) {
    printf(" Peigs: ci_size_  node %d : Mapping problem\n", *me);
    indx = -3;
    xerbla_("ci_siz\n", &indx);
  }
  
  jndx = -3;
  
  ptr = map;
  for ( indx = 0; indx < *n; indx++ ) {
    if ( (ptr++) == NULL )
      xerbla_("ci_siz\n", &jndx);
  }
  
  jndx = 0;
  ptr = map;
  kndx = *n;
  iam = *me;
  for ( indx = 0; indx < *n; indx++ ) {
    if ( *(ptr++) == iam )
      jndx += kndx - indx;
  }
  
  return(jndx);
}

Integer fil_mapvec_( me, n, map, mapvec)
     Integer *me, *n, *map, *mapvec;
{
  /*
    from the map array construct a shorter
    array with information about the vectors stored on this processor
    
    mapvec[i] = j
    i = the i-th vector stored on this processor
    j = the real column/row indx of the i-th vector
    
    returns the number of vectors on this processor given
    by map and also the mapvec list
    
    argument:

    me         = (input/integer) node id
    
    n           = (input/integer) dimension of the matrix
    
    map       = (input/integer array) distribution of the columns of the matrix
    				length n
    
    mapvec  = (output/integer array) length n
    
    */
  
  Integer i, k, *iptr, idummy;
  
  extern Integer mxnprc_();
  extern Integer mxmynd_();
  
  extern void xerbla_();

  if ( me == NULL ) {
    i = -1;
    printf(" Peigs: fil_mapvec_  node %d : first argument is NULL \n", mxmynd_());
    printf("me = %d fil_mapvec_  \n", *me);
  }
  
  if (( *me > mxnprc_() ) || ( *me < 0 )) {
    i = -1;
    printf(" Peigs: fil_mapvec_  node %d : first argument %d out of bounds \n", mxmynd_(), *me);
    printf("fil_mapvec_ %d \n", i);
  }
  
  if ( n == NULL ) {
    i = -2;
    printf(" Peigs: fil_mapvec_  node %d : 2nd argument is NULL i\n", *me);
    printf("fil_mapvec_ %d \n", i);
  }
  else if ( *n < 1 ) {
    i = -2;
    printf(" Peigs: fil_mapvec_  node %d : Mapping problem second argument is invalid \n", *me);
    printf("fil_mapvec_ %d \n", i);
  }
  
  if (map == NULL) {
    printf(" Peigs: fil_mapvec_  node %d : Mapping problem\n", *me);
    i = -3;
    printf("fil_mapvec_ %d \n", i);
  }
  
  iptr = map;
  for ( i = 0 ; i < *n ; i++ )
    if ( (iptr++) == NULL ) {
      i = -3;
      printf(" Peigs: fil_mapvec_  node %d : 3rd argument error. \n", *me);
      printf("fil_mapvec_ %d \n", i);
    }
  
  if (mapvec == NULL) {
    printf(" Peigs: fil_mapvec_  node %d : Mapping problem\n", *me);
    printf(" Peigs: fil_mapvec_  node %d : 4th argument error. \n", *me);
    i = -4;
    printf("fil_mapvec_ %d \n", i);
  }
  
  
  k = 0;
  iptr = map;
  idummy = *me;
  for ( i = 0; i < *n; i++ ) {
    if ( *(iptr++) == idummy ) {
      if ( mapvec + k == NULL ) {
        i = -4;
	printf(" Peigs: fil_mapvec_  node %d : 4th argument error. \n", *me);
	printf("fil_mapvec_ %d \n", i);
      }
      mapvec[k] = i;
      k++;
    }
  }
  
  return(k);
}


