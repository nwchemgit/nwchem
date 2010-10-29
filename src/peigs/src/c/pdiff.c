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
#include <stdlib.h>

#include "globalp.c.h"

void pdiff( n, data, proclist, nprocs, iwork, msg1, msg2, info )
     char    *data, *msg1, *msg2;
     Integer *n, *proclist, *nprocs, *iwork, *info;
{
/*
    This routine performs an element by element comparision of a
    length n array of characters, data, on a processor with its
    neighbor as described by a list of processors proclist[0:nprocs-1].
    
    n ........... number of bytes of data in array data

    *data ....... pointer to char (char *)
                  character array with n elements which should be the
                  same on all processors in proclist

    *nprocs ...... number of processor ids in proclist

    *proclist ... pointer to Integer (Integer *)
                  a list of 'nprocs' processor ids on which 'data' should be
                  the same.

                  There must be no repeated processor ids in proclist
   
    iwork ....... pointer to Integer (Integer *)
                  workspace of length at least n bytes, i.e., the same
                  size as data
 
    msg1 ........ pointer to char (char *)
                  General error message, e.g. "Error in pdspevx"

    msg2 ........ pointer to char (char *)
                  Specific error message indicating what data represents
                  so can tell user what data differs on different
                  processors

    *info = 0:  'data' on this processor is the same as 'data' on
                one of this processor's neighbors.

          = 1:  'data' on this processor is different from 'data' on
                one of this processor's neighbors.  *info = 1 on one
                processor does not mean that *info != 0 on all processors
                in proclist.  To do this you need to do a global operation.
                This is typically by pgexit().

                The processors on which *info = 1 print an error message
                to stderr
 */
  
  static Integer TYPE = 10;

  Integer isize, me, last_proc, next_proc, indx;
  
  extern Integer indxL ();
  extern Integer mxwrit_(), mxread_(), mxmynd_();
  

  *info = 0;

  me = mxmynd_();

  indx = indxL ( me, *nprocs, proclist);
  
  if ( indx < 0 ) {
    fprintf( stderr, " Error in pdiff.  me = %d not in proclist. \n", me );
    exit(-1);
  }

  if ( *nprocs < 2 )
     return;
    
  last_proc = (indx + *nprocs - 1) % *nprocs;
  last_proc = proclist[last_proc];
  next_proc = (indx + 1) % *nprocs;
  next_proc = proclist[next_proc];
    
  isize = *n;

  if ( ( (indx + 2) % 2) == 0 ) {
    mxwrit_( (Integer *)data,  &isize, &next_proc, &TYPE );
    mxread_( iwork, &isize, &last_proc, &TYPE );
  }
  else {
    mxread_( iwork, &isize, &last_proc, &TYPE );
    mxwrit_( (Integer *)data,  &isize, &next_proc, &TYPE );
  }
    
  indx = memcmp( data, iwork, *n );
  indx = abs(indx);

  if( indx != 0 ) {
    fprintf( stderr, " %s %s is different on processors %d and %d \n",
             msg1, msg2, me, last_proc );

    *info = 1;
  }

  return;
}

