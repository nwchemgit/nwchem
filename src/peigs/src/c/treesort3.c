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
  peigs sort routine: tree sort on an integer array 
  */ 


#include <stdio.h>
#include <memory.h>
#include <stdlib.h>

#include "globalp.c.h"

static Integer peigs_malloc_adr;

void gtreeinsert(p, list, iscrat)
     Integer p;
     Integer *list, *iscrat;
{
	Integer i;

  if ( p!= -1 ){
  	gtreeinsert(iscrat[p+2], list, iscrat);
	i = iscrat[p+1];
  	list[i] = iscrat[p];
  	gtreeinsert(iscrat[p+3], list, iscrat);
  }
  return;
}

Integer gtree(p, w, count, iscrat)
     Integer p;
     Integer w, *count, *iscrat;
{
  extern Integer mxmynd_ ();
  Integer gtree();

if ( p == -1 ) {
	p = peigs_malloc_adr;
	peigs_malloc_adr += 4;
	iscrat[p] = w;
  	iscrat[p+1] =  *count;
	iscrat[p+2] = -1;
	iscrat[p+3] = -1;
	(*count)++;
}
  else {
    if ( w < iscrat[p] )
      iscrat[p+2] = gtree( iscrat[p+2], w, count, iscrat);
    else if ( w > iscrat[p] )
      iscrat[p+3] = gtree( iscrat[p+3], w, count, iscrat);
  }
  
  return(p);
}

Integer reduce_list3( num_list, list, scratch )
     Integer num_list, *list, *scratch;
{
  /*
    this is just a quick hack
    
    ---should use some hash or avl-tree for better long term solution
    
    this utility removes duplicate elements in the list
    the content of list and returns the "shorter list" in scratch
    
    input
    
    num_list = number of element in list
    
    output = the number of distinct elements in scratch
    preserving the order in which they appear
    
    */
  
  Integer *junk, k, num_indx;
  Integer root, gtree(), *iscrat;

  peigs_malloc_adr = 0;
  if ( num_list == 0 )
    return(0);

  root = -1;
  junk = scratch;

/*
  iscrat = junk + num_list;
*/

  iscrat = (Integer *) malloc(5*num_list*sizeof(Integer));
  if ( iscrat == NULL ) {
        fprintf(stderr, " Not enough memory for junk; malloc fails \n");
        return(-1);
        }
  
  for ( k = 0; k < 5*num_list; k++ )
    iscrat[k] = -1;
  
  num_indx = 0;
  for ( k = 0; k < num_list; k++ )
    root = gtree(root, list[k], &num_indx, iscrat);
  
  for ( k = 0; k < num_indx; k++ )
    scratch[k] = -1;
  
  gtreeinsert(root, scratch, iscrat);

  free(iscrat);
  return(num_indx);
}

Integer reduce_list33( num_list, list, scratch, node )
     Integer num_list, *list, *scratch, node;
{
  /*
    
    this is just a quick hack; replace reduce_list
    
    ---should use some hash or avl-tree for better long term solution
    
    this utility removes duplicate elements in the list
    the content of list and returns the "shorter list" in scratch
    
    input
    
    num_list = number of element in list
    
    output = the number of distinct elements in scratch
    preserving the order in which they appear

    */

  Integer *junk, k, num_indx;
  Integer root, gtree();

  if ( num_list == 0 )
    return(0);

  root = -1;
  num_indx = 0;
/*
  junk = scratch + num_list;
*/
  junk = (Integer *) malloc(5*num_list*sizeof(Integer));
  if ( junk == NULL ) {
	fprintf(stderr, " Not enough memory for junk; malloc fails \n");
	return(-1);
	}

  for ( k = 0; k < 5*num_list; k++ )
    junk[k] = -1;

  root = gtree(root, node, &num_indx, junk);

  for ( k = 0; k < num_list; k++ )
    root = gtree(root, list[k], &num_indx, junk);
  
  for ( k = 0; k < num_indx; k++ )
    scratch[k] = -1;
  
  gtreeinsert(root, scratch, junk);
  
  free(junk);
  
  return(num_indx);
}




Integer reduce_list2( num_list, list, scratch )
     Integer num_list, *list, *scratch;
{
  /*
    this is just a quick hack
    
    ---should use some hash or avl-tree for better long term solution
    
    this utility removes duplicate elements in the list
    the content of list and returns the "shorter list" in scratch
    
    input
    
    num_list = number of element in list
    
    output = the number of distinct elements in scratch
    preserving the order in which they appear
    
    */
  
  Integer k, i, last_num;
  Integer gtree();

  if ( num_list == 0 )
    return(0);
  
  for ( i = 0; i < num_list; i++ )
    scratch[i] = list[i];
  k = qqsort(scratch, 0, num_list-1);
  
  k = 0;
  last_num = scratch[0];
  
  /*
    scratch is sorted 
    */
  
  for ( i = 1; i < num_list; i++ ) {
    if ( last_num != scratch[i]){
      k++;
      last_num = scratch[i];
      scratch[k] = scratch[i];
    }
  }
  
  k++;  
  return(k);
}


Integer reduce_list4( num_list, list, list2, iwork )
     Integer num_list, *list, *list2, *iwork;
{
  /*
    
    this utility removes duplicate elements in the list
    the content of list and returns the "shorter list" in list2
    iwork must be of length mxnprc_().
    
    input
    
    num_list = number of element in list
    list     = list of processor ids to be reduced.
    
    !!! MUST have 0 < list[i] < mxnprc_().  This is not checked !!!
    
    output = the number of distinct elements in iwork  
    preserving the order in which they appear
    
    */
  
  Integer naproc, i, nproc;
  extern Integer mxnprc_();

  if ( num_list == 0 )
    return(0);
  
  naproc = mxnprc_();
  
  for ( i = 0; i < naproc; i++ )
    iwork[i] = 0;

  nproc = 0;
  for ( i = 0; i < num_list; i++ )
    if( iwork[ list[ i ] ] == 0 ) {
      list2[ nproc ] = list[ i ];
      iwork[ list[ i ] ] = 1;
      nproc++;
    }
  
  return(nproc);
}

Integer reduce_list22( num_list, list, scratch, node )
     Integer num_list, *list, *scratch, node;
{
  /*
     scratch is length n+1
     
     */

  Integer k, i, last_num;
  Integer gtree();
  
  if ( num_list == 0 )
    return(0);
  
  for ( i = 0; i < num_list; i++ )
    scratch[i] = list[i];
  scratch[num_list] = node;
  k = qqsort(scratch, 0, num_list);
  
  k = 0;
  last_num = scratch[0];
  
  /*
     scratch is sorted 
     */
  
  for ( i = 1; i < num_list+1; i++ ) {
    if ( last_num != scratch[i] ){
      k++;
      last_num = scratch[i];
      scratch[k] = scratch[i];
    }
  }
  
  k++;  
  return(k);
}




