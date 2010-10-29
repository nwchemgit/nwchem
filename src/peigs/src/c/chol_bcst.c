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
/* *********************************************
   
   PeIGS communication routine
   
   pipeline broadcast;
   
   name because we first used
   it in the choleski factorization 
   
   pipeline broadcast from root
   all duplicate processor ids
   are removed from the list of processors call map.
   
   chol_pipe_bcast(buf, len, type, root, k_indx, map, scratch)
   
   status: not protected from user input errors
   
   */


#include <stdio.h>
#include <memory.h>

#include "globalp.c.h"

void chol_pipe_bcast(buf, len, type, root, k_indx, map, scratch)
     char buf[];
     Integer len;
     Integer type;
     Integer root, k_indx;
     Integer map[], scratch[];
     
/*

in this version of broadcast, we are broadcasting to a list of
processors, proclist[0:k_indx-1] from root
all duplicate node ids are removed so that we have a
minimal set of processors doing the pipelining

buf         = (input/char) this is the buffer that is to
                 be broadcasted ( should move to untyped)

len          = (input/integer) this is the size, in bytes,
                  of the buffer

type        = (input/integer) type for the message type

root         = (input/integer) the root id

k_indx   = (input/integer) length of the list in map

map        = (input/integer array) is a ptr to the
beginning of a
list of processor.  Typically this is the
locations of the columns/rows
		     of the matrices
		     
scratch   = (integer scratch array) length(k_indx + 1)

*/
{
  Integer *i_ptr;
  Integer num_indx, me;
  Integer *junk, k, isize;
  
  extern Integer in_list();
  extern Integer qqsort();
  extern void qshellsort_ ();
  
  extern Integer mxmync_ ();
  extern Integer mxmynd_ ();
  extern Integer mxwrit_ ();
  extern Integer mxread_ ();
  extern Integer reduce_list22();
  
  /*
    reduce the list of processors to the minimum list holding 
    minimal number of processors
    */
  
  me = mxmynd_ ();
  
  if ( k_indx < 1 )
    return;

  junk = scratch;  

  junk[0] = root;
  num_indx = 1;
  for ( k = 0; k < k_indx; k++ ){
  i_ptr = &map[k];
  if ( ( map[k] != root ) && ( in_list( i_ptr, junk, &num_indx) == 0 )){
    junk[num_indx] = map[k];
    num_indx++;
  }
}
  
/*  
  num_indx = reduce_list22(k_indx, map, junk, root);
*/
  
  /*
    junk[0:num_indx] now contains the list of processors in order
    without root and duplication
    in sequence of the method they appear according to map[*]
    */
  
  if ( num_indx == 1 ) 
    return;  /* root is the only node */
  
  /*
    quick sort the processor list;
    this may improve the communication
    on some architecture; should
    incorporate a variant of quick sort
    and a reducelist in one
    
    qqsort( junk, 1, num_indx-1);
    */
  
  isize = num_indx-2;
  
#ifdef MESH
  gshellsort_ ( &isize, &junk[1]);
#endif
  
  if ( me == root ) {
    isize = mxwrit_ ( buf, &len, &junk[1], &type);
  }
  
  for ( k = 1; k < num_indx; k++ ){
    if ( me == junk[k] ) {
      isize = mxread_ ( buf, &len, &junk[k-1] , &type);
      if ( k < ( num_indx - 1 ) ){
	isize = mxwrit_ ( buf, &len, &junk[k+1] , &type);
      }
    }
  }
  return;
}


