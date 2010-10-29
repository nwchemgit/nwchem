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
  PeIGS internal error checking routine.

  void mapdif_( n, mapA, mapB, iscratch, ndiff )
  Integer            *n, *mapA, *mapB, *iscratch, *ndiff;

  Compares the set of processors defined by mapA
  and the set of processors defined by mapB.
  
  
  */

#include <stdio.h>
#include <string.h>

#include "globalp.c.h"

void mapdif_( n, mapA, mapB, iscratch, ndiff )
     Integer            *n, *mapA, *mapB, *iscratch, *ndiff;
     
     /*
      *  Routine to determine if mapA and mapB contain the same set of processor ids.
      *
      *  On Entry:
      *    n .......... Number of elements in mapA and mapB
      *    mapA ....... A list of n processor ids.  All ids must be 0 <= id < mxnprc_().
      *    mapB ....... B list of n processor ids.  All ids must be 0 <= id < mxnprc_().
      *    iscratch  .. Workspace of length ( 2 * n )
      *
      *  On Exit:
      *    ndiff ...... = 0 then the set of unique processor ids in mapA and mapB are the 
      *                     same (though the ids may not be in the same order or repeated
      *                     the same number of times.)
      *                 > 0 then mapA has processor ids which are not included in mapB and/or
      *                          mapB has processor ids which are not included in mapA.
      *                         
      */
     
{
  Integer             nproc;
  Integer             *ptrA, *ptrB;
  Integer             *iscrat;
  Integer nprocA, nprocB, status;
  Integer me;
  
  extern Integer      mxnprc_();
  extern Integer      mxmynd_();
  extern Integer      reduce_list2();
  extern Integer      qqsort();
  
  nproc = mxnprc_();
  me = mxmynd_ ();
  
  iscrat = iscratch;
  ptrA = iscrat;
  
  nprocA = reduce_list2( *n, mapA, ptrA );
  iscrat += nprocA;
  ptrB = iscrat;
  nprocB = reduce_list2( *n, mapB, ptrB );
  iscrat += nprocB;
  
  if (( nproc >= nprocA ) && (nproc >= nprocB )) {
    if ( nprocB == nprocA ) {
      qqsort(ptrA, 0, nprocA-1);
      qqsort(ptrB, 0, nprocB-1);
      status = nprocB * sizeof(Integer);
      *ndiff = memcmp( ptrA, ptrB, status);
      return;
    }
    else {
      *ndiff = -10;
      return;
    }
  }
  else {
    *ndiff = -100;
    return;
  }
}
