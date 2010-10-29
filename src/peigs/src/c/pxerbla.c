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

PeIGS internal error utility: pxerbla2_

global list check of some a given list

not intended for separate call or protected from
user input errors

*/

#include <stdio.h>
#include <memory.h>
#include <math.h>
#include <stdlib.h>

#include "globalp.c.h"

void pxerbla2_( n, array, procmap, len, iwork, info )
     char *array;
     Integer *n, *procmap, *len, *iwork, *info;
{
  /*
    This routine performs an
    element by element comparision of a length n array of characters on
    a processor with its neighbor as described by a list of processors
    map[0:len-1].
    
    It is assumed that values in the array map[0:len-1] form
    the same set of processors (User Beware); for this routine, the
    ordering in map[0:len-1] is not important.  if the array[0:n-1] on
    each of the processors in map[0:len-1] does not match pxerbal1_
    returns

    *info > 0:  array on this processor is different from "array" on
                one of this processor's neighbor.  *info != 0 on one
                processor does not mean that *info != 0 on all processors
                in procmap.  To do this you need to do a global operation.
                This is typically by g_exit_().
    
    *info = -50: this processor is not in procmap
    *info = -1:  number of distinct processor ids in procmap > # of 
                 allocated processors
    
    arguments:
    
    n = number of bytes in a character array[*]

   WORKSPACE
    let nproc = Number of unique processor ids in procmap, i.e.,
                nprocs = reduce_list( *len, procmap, proclist).

    Then:
 
    iwork   = scratch array of length ( nproc + room for *n char variables)

    */
  
  static Integer TYPE = 10;
  Integer isize, nprocs, me, me_indx, maxprocs;
  Integer last_proc, next_proc, indx;
  
  Integer *iscrat, *proclist, *map_in;
  
  extern Integer reduce_list2();
  extern Integer indxL ();
  extern Integer mxwrit_(), mxread_();
  extern Integer qqsort();
  extern void gi_sum();
  extern void xerbla_();
  extern Integer mxcmp();
  extern Integer      menode_();
  extern Integer mxmynd_();
  extern Integer mxnprc_();
  

  *info = 0;

  me = mxmynd_();

#ifdef DEBUG2
  fprintf(stderr, " in pxerbla me = %d \n", me);
#endif

  maxprocs = mxnprc_();       /* the maximum number of processors allocated */
  
  iscrat = iwork;
  proclist = iscrat;
  nprocs = reduce_list2( *len, procmap, proclist);
  iscrat += nprocs;

  qqsort( proclist, 0, nprocs-1);
  
  if ( nprocs > maxprocs ) {
    fprintf(stderr, "PXERBLA: Node %d Error: Number of processors in Proc List exceeds number allocated \n", me);
    *info = -1;
    return;
  }
  
  /*
    the number of distinct processors
    */
  
  map_in = iscrat;
  indx = me_indx = indxL ( me, nprocs, proclist);
  
  /*
    i am actually in the list participating in this check
    */
  
  if ( indx != -1) {
    if ( nprocs == 1 )
      return;
    
    last_proc = (me_indx + nprocs - 1) % nprocs;
    last_proc = proclist[last_proc];
    next_proc = (me_indx + 1) % nprocs;
    next_proc = proclist[next_proc];
    /*
      isize = *n * sizeof(Integer);
      */
    
    isize = *n;
    if ( (me_indx % 2) == 0 ) {
      mxwrit_( array, &isize, &next_proc, &TYPE );
      mxread_( map_in, &isize, &last_proc, &TYPE );
    }
    else {
      mxread_( map_in, &isize, &last_proc, &TYPE );
      mxwrit_( array, &isize, &next_proc, &TYPE );
    }
    
    isize = *n;
    indx = memcmp( array, map_in, isize );
    indx = abs(indx);
    *info = indx;
    return;
  }
  else
    *info = -50;  /* not even in the list of processors */

#ifdef DEBUG2
  fprintf(stderr, " out pxerbla2 me = %d \n", me);
#endif

  return;
}

