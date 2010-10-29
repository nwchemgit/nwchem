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

void mdif2b_( nA, mapA, nB, mapB, nC, mapC, proclist, ndiff )

  Integer      *nA, *mapA, *nB, *mapB, *nC, *mapC, *proclist, *ndiff;
     
/*
 *  Returns a list of any processor ids in mapA but not in mapB and/or mapC.
 *  i.e., proclist = {mapA} - ( {mapB}) U {mapC} )
 *
 *  ON ENTRY
 *  --------
 *
 *    nA ......... Number of elements in mapA
 *
 *    mapA ....... A list of nA processor ids.  All ids must be 0 <= id < mxnprc_().
 *
 *    nB ......... Number of elements in mapB
 *
 *    mapB ....... A list of nB processor ids.  All ids must be 0 <= id < mxnprc_().
 *
 *    nC ......... Number of elements in mapC
 *
 *    mapC ....... A list of nC processor ids.  All ids must be 0 <= id < mxnprc_().
 *
 *    proclist  .. Vector of length ( mxnprc_() )
 *
 *  ON EXIT
 *  -------
 *
 *    ndiff ...... = 0 then all processor ids appearing in mapA also appear
 *                          mapB and/or mapC.
 *                 > 0 then mapA has processor ids which do not appear in
 *                          mapB and/or mapC.
 *
 *    proclist ... If ndiff = 0: then junk
 *                 Otherwise:    The first ndiff elements is a list of the processor
 *                               ids (without any repeats)
 *                               which appear in mapA, but not mapB and/or mapC, sorted
 *                               in increasing order.
 *                         
 */
     
{
  Integer             nproc, indx, ndiff1;
  
  extern Integer      mxnprc_();
  


  nproc = mxnprc_();
    
  for (indx = 0; indx < nproc; indx++ )
     proclist[indx] = 0;
    
  for (indx = 0; indx < *nA; indx++ )
     proclist[ mapA[ indx ] ] = 1;
    

  for (indx = 0; indx < *nB; indx++ )
     proclist[ mapB[ indx ] ] = 0;
    
  for (indx = 0; indx < *nC; indx++ )
     proclist[ mapC[ indx ] ] = 0;
    

  ndiff1 = 0;
  
  for (indx = 0; indx < nproc; indx++)
    if( proclist[ indx ] > 0){
      proclist[ ndiff1 ] = indx;
      ndiff1++;
    }
  
  for (indx = ndiff1; indx < nproc; indx++)
    proclist[ indx ] = -1;
  
  *ndiff = ndiff1;
  
  return;
}
