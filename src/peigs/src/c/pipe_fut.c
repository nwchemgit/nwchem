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
/* *******************************************
   
   
   PeIGS communication routine
   
   void pipe_bcst_fut_col(n, buf, len, type, root, c_indx, map, scratch)
   
   pipeline communication to a list of processors with index
   bigger than c_indx

   *
   *
   */


   
#include <stdio.h>
#include <math.h>

#include "globalp.c.h"

void pipe_bcst_fut_col(n, buf, len, type, root, c_indx, map, scratch)
     Integer n;
     char *buf;
     Integer len;
     Integer type;
     Integer root, c_indx;
     Integer *map, *scratch;
     
     /* 
       
       in this version of broadcast, we are broadcasting to a list of
       processors holding columns with column numbers greater than c_indx
       
       duplications in the processor list is removed
       
       c_indx  = current index; the array map is assumed to have size at least
       
       
       map     = is a ptr to the beginning of a
       list of locations of the columns/rows of the matrices
       
       root    = the node that started it all
       
       scratch = scratch integer memory of length k_indx
       
       ;;;should probably sort the list for maximum efficiency on the delta and hypercubes
       
       */
{
  Integer num_indx, me;
  Integer *junk, *ptr, k, isize;
  extern Integer mxmynd_();
  extern Integer mxwrit_();
  extern Integer mxread_();
  extern Integer in_list();

/*
  reduce the list of processors to a minimum list holding previous column numbers
*/

  me = mxmynd_ ();
  
  if ( c_indx < 0 )  /* the 0-th column has no previous columns */
    return;
  
  junk = scratch;
  *junk = root;
  num_indx = 1;
  
  /*
    reduce the list to only the shortest processor list needed
    */
  
  for ( k = c_indx+1; k < n; k++ ){
      if ( ( map[k] != root ) && ( in_list( &map[k], junk, &num_indx) == 0 )){
	*(junk + num_indx) = map[k];
	num_indx++;
      }
    }
  
  /*
    num_indx holds the number of processors holding future columns
    */
  
  /*
    junk[1:num_indx] now contains the list of processors without root and duplication
    in sequence of the method they appear according to map[*]
    */
  
  if ( num_indx == 1 ) 
    return;
  
  /* root is the only node */
  
  ptr = &junk[1];

  if ( me == root ) {
    isize = mxwrit_( buf, &len, ptr, &type);
  }
  
  for ( k = 1; k < num_indx; k++ ){
    if ( me == *ptr ){
      isize = mxread_ ( buf, &len, ptr-1 , &type);
      if ( k < ( num_indx - 1 ) )
	isize = mxwrit_( buf, &len, ptr+1 , &type);
    }
    ptr++;
  }
  return;
}
