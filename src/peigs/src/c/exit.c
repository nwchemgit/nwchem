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
/* *****************************************
   
   PeIGS internal utility
   
   g_exit_
   
   global list error exit check
   
   */


#include <stdio.h>
#include <memory.h>
#include <math.h>
#include <stdlib.h>

#include "globalp.c.h"


void g_exit2_( n, array, procmap, len, iwork )
     char *array;
     Integer *n, *procmap, *len, *iwork;
{
  /*

      An old, out of date version of g_exit_.

      This routine should never be called, use g_exit_ instead.  To use
      g_exit2_ in place og g_exit2_ you just need change the name
      of g_exit2 to g_exit, and add a DoublePrecision workspace containing
      at least bufsiz bytes (see cmbbrf.h) to the end of the argument list.
      

      If this routine is called then it prints an error
      message and aborts program execution.
    
    */
  
  Integer me;
  extern void mxpend_ ();
  extern Integer mxmynd_();

  me = mxmynd_ ();
  
  fprintf(stderr, "G_EXIT2: Node %d Error.  A routine called g_exit2_, \n", me);
  fprintf(stderr, "Message from calling routine: %s  node id = %d \n", array, me);

  mxpend_ ();

  exit(-1);

  return;
}


void g_exit_( n, array, procmap, len, iwork, work )
     char *array;
     Integer *n, *procmap, *len, *iwork;
     DoublePrecision *work;  /* workspace containing at least bufsiz bytes (see cmbbrf.h) */
{
  /*
    This routine performs a global check on an integer n: if n is negative
    on any processor in procmap the routine exits with a message from array.
    n should be less than or equal to 0 on all processors in procmap.
    
    a global combine is performed on n.
    if n is less than 0 upon return, a exit is called.
    
    It is assumed that values in the array map[0:len-1] form the same set
    of processors (User Beware); for this routine, the ordering in map[0:len-1]
    is not important.
    
    if n is not 0 after a global combine then all processors in procmap
    set *n = -51, print the error message in "array", plus a statement
    indicating that info = -51 and exits.
    
    n       = integer
    array   = character string for messages
    procmap = array of processors on which to check n
    len     = length of the array procmap

   WORKSPACE
    let nproc = Number of unique processor ids in procmap, i.e.,
                nprocs = reduce_list( *len, procmap, proclist).

    Then:
 
    iwork   = scratch array of length ( *len )

    
    */
  
  static Integer TYPE = 10;
  Integer nprocs, me, maxprocs;
  Integer *iscrat, *proclist;
  
  extern Integer reduce_list2();
  extern Integer indxL ();
  extern Integer mxwrit_ (), mxread_ ();
  extern Integer qqsort();
  extern void gi_sum();
  extern void xerbl2_ ();
  extern Integer mxcmp();
  extern void mxpend_ ();
  extern Integer mxmynd_();
  extern Integer mxnprc_();

  me = mxmynd_ ();
  maxprocs = mxnprc_ ();       /* the maximum number of processors allocated */
  
  iscrat = iwork;
  proclist = iscrat;
  nprocs = reduce_list2( *len, procmap, proclist);
  iscrat += nprocs;
  qqsort( proclist, 0, nprocs-1);
  
  if ( nprocs > maxprocs ) {
    fprintf(stderr, "G_EXIT: Node %d Error: Number of processors in Proc List exceeds number allocated \n", me);
    xerbl2_ ( );
    return;
  }
  
  /*
    i am actually in the list participating in this check
    */
  
  gi_sum( n, 1, TYPE, proclist[0], nprocs, proclist, work);
  
  if ( *n < 0 ) {
    *n = -51;
    fprintf(stderr, " %s  My node id = %d info = %d (g_exit_) \n", array, me, *n);
    /*
      xerbl2_ ( );
      */
    mxpend_ ();
    exit(-1);
  }
  return;
}


void gi_sum(buf, items, msgtype, root, snumprocs, plist, work)
     /*
       this is a integer global sum on buf
       */
     Integer *buf;
     Integer items;
     Integer msgtype;
     Integer root;
     Integer *plist, snumprocs;
     DoublePrecision *work;  /* workspace containing at least bufsiz bytes (see cmbbrf.h) */
{
  Integer isize;
  extern Integer sumiv_();
  extern Integer mxcombv1_();
  
  isize = sizeof(Integer);

  mxcombv1_ ( (char *) buf, sumiv_ , &isize, &items, &snumprocs, plist, &msgtype, (char *)work);

  return;
}

  

