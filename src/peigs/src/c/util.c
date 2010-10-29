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
/* 
  this contains a number of the utility programs used by peigs to get
  information about lists:

  count_list  : count number of occurrance of certain integer in a list

  indxL       : get the index of some element in a list

  indxlf_     :  fortran callable version of indxl

  indaint     : this routine returns the first index indaint of
                the array map such that map[indaint] >= k


  fil_int_lst  : fill an integer array with integers

  fil_dbl_lst  : fill an array doubles with a DoublePrecision

  in_list      : is some integer in a list

  find_proc_store:  from a list of natural numbers find the first index with non-zero

  find_large_store:  from a list of natural numbers find the first index with largest non-zero number

  mem_cpy        : integer copy


*/


#include <stdio.h>
#include <math.h>

#include "globalp.c.h"

#define max(a,b) ((a) > (b) ? (a) : (b))
#define min(a,b) ((a) < (b) ? (a) : (b))

Integer count_list ( me, list, size )
     Integer me, list[], *size;
{
  /*
    
    count the number of instance of "me" in a "list"
    
    */
  
   Integer i, j;
   Integer *ptr;
	Integer isize  = *size;
  
  if ( isize <= 0 )
    return(0);
  
  ptr = list;
  j = 0;  
  for ( i = 0; i < isize; i++ ) {
    if ( *(ptr++) == me  ) 
      j++;
  }
  return(j);
}


Integer indxL( k, nvecs, map)
     Integer k, nvecs, map[];
     /*
       this routine computes the index of the real column number
       
       e.g. say real column number is k
       map[indx] = k
       
       this routine returns the first indx that it encounters
       
       */
{
  Integer i;
  Integer *ptr;
  
  ptr = map;
  for ( i = 0; i < nvecs; i++ ){
    if ( ptr[i] == k )
      return(i);
  }
  return(-1);  /* failed */
}


Integer indxlf_ ( k, nvecs, map)
     Integer *k, *nvecs, map[];
     /*
       this routine computes the index of the real column number
       
       e.g. say real column number is k
       map[indx] = k
       
       this routine returns indx
       
       */
{
  Integer i;
  Integer *ptr;
	Integer nnvecs = *nvecs;
	Integer kk = *k;
  
  ptr = map;

  for ( i = 0; i < nnvecs; i++ )
    {
      if ( *(ptr++) == kk )
	return(i);
    }
  return(-1);  /* failed */
}

Integer indaint( k, nvecs, map)
     Integer k, nvecs, map[];
     /*
       this routine returns the first index indaint of the array map
       such that map[indaint] >= k
       
       */
{
  Integer i;
  Integer *ptr;
  
  ptr = map;
  for ( i = 0; i < nvecs; i++ )
    {
      if ( *(ptr++) >= k )
	return(i);
    }
  return(-1);  /* failed */
}

void fil_int_lst ( n, list, item)
     Integer n, item;
     Integer list[];
{
  /*
    fill an integer array list with item
    */
  
  Integer i, w;
  Integer *j;
  
  j = list;
  w = item;
  for ( i = 0; i < n; i++ ) {
    *(j++) = w;
  }
  return;
}


void fil_dbl_lst ( n, list, item)
     Integer n;
     DoublePrecision list[], item;
     
{
  /*
    fill an DoublePrecision precision array list with item; is there an assembly version
    of this?
    */
  
   Integer i;
   DoublePrecision w;
  
  w = item;
  for ( i = 0; i < n; i++ ) {
    list[i] = w;
  }
  return;
}

Integer in_list( item, list, list_len)
     Integer *item, list[], *list_len;
{
  /*
    this routine looks through the *list of length *list_len
    for *item
    
    returns 1 if found
    0 if not
    */
  
  Integer i;
  Integer tmp;
  Integer len;
  
  tmp = *item;
  len = *list_len;
  for ( i = 0; i < len; i++ ) {
    if( tmp == list[i] )
      return(1);
  }
  return(0);
}

Integer find_proc_store( indx, nprocs, size, sizelist)
     Integer indx, nprocs,  size, sizelist[];
{
  /*
    returns the first integer $i in sizelist such that
    sizelist[ $i ] > size
    */
  
  Integer i, j;
  
  for ( i = indx ; i < nprocs + indx; i++ ) {
    j = i % nprocs;
    if ( sizelist[j] > size )
      return(i);
  }
  return(-1);
}
	

Integer find_large_store( indx, nprocs, size, sizelist)
     Integer indx, nprocs,  *size, sizelist[];
{
  /*
    returns the last integer $i in sizelist
    with the largest amount which can accomodate size
    
    if no such integer $i existst then return the indx
    with the largest size less than size

    indx =
    nprocs = number of processors in the list
    size = number of vectors
    sizelist = list of vector memory available on participating processors
    
    */
  
  Integer i, j, k;
  Integer lsize, id, l,m;

  lsize = *size;
  
  j = -1;
  k = -1;
  l = -1;
  m = -1;

  for ( i = indx ; i < (nprocs + indx); i++ ) {
    id = (i % nprocs);
    if ( sizelist[id] >= lsize ) {
      l = max( l, sizelist[id] );
      if ( l == sizelist[id])
	m = id;
    }
    else
      {
	k = max(k, sizelist[id]);
	if ( k == sizelist[id] )
	  j = id;
      }
  }
  /*
    l = largest size bigger than lsize
    k = largest size overall
    */

  if ( l == -1 ) {
    *size = min(k, *size);
    return(j);
  }
  else {
    *size = min(l, *size);
    return(m);
  }
}

void mem_cpy ( list1, list2, n)
     Integer list1[], list2[], n;
     /*
       
       copies list1[0:n-1] to list2[0:n-1]
       
       should compare with memcpy
       
       */
{
  Integer i;
  Integer *ptr1, *ptr2;
  
  
  ptr1 = list1;
  ptr2 = list2;
  
  for ( i = 0; i < n; i++ ){
    *(ptr2++) = *(ptr1++);
  }
  
  return;
}
