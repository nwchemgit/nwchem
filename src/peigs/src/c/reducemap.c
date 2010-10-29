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

void reduce_maps( nA, mapA, nB, mapB, nC, mapC, nproc, proclist )

  Integer      nA, nB, nC;
  Integer      mapA[], mapB[], mapC[], *nproc, proclist[];
     
/*
 *  Combine the processor ids in mapA, mapB, and mapC into a single sorted
 *  listed without duplication.  May reduce/combine just one or two maps
 *  by setting two, or one, respectively of nA, nB, and nC to zero.
 *
 *    nA ......... (Input) Integer
 *                 Number of elements in mapA, >= 0.
 *
 *    mapA ....... (Input) pointer to Integer array of length (n) = (Integer *)
 *                 A list of nA processor ids.  All ids must be 0 <= id < mxnprc_().
 *
 *    nB ......... (Input) Integer
 *                 Number of elements in mapB >= 0.
 *
 *    mapB ....... (Input) pointer to Integer array of length (n) = (Integer *)
 *                 B list of nB processor ids.  All ids must be 0 <= id < mxnprc_().
 *
 *    nC ......... (Input) Integer
 *                 Number of elements in mapC >= 0.
 *
 *    mapC ....... (Input) pointer to Integer array of length (n) = (Integer *)
 *                 C list of nZ processor ids.  All ids must be 0 <= id < mxnprc_().
 *
 *    nproc ...... (Output) pointer to Integer = (Integer *)
 *                 Number of distinct processor ids in combined mapA, mapB and mapC.
 *
 *    proclist ... (Output) pointer to Integer array of length (mxnprc_() >0 nproc) = (Integer *)
 *                 Vector containing the nproc unique processor
 *                 ids in combined mapA, mapB and mapC in sorted order.
 *
 *                 Note that proclist[nproc:mxnprc_()-1] is used as workspace
 *                 
 *
 */
     
{
  Integer             indx, naproc, mproc;
  
  extern Integer      mxnprc_();
  

  naproc = mxnprc_();
    
  for (indx = 0; indx < naproc; indx++ )
     proclist[indx] = 0;
    
  for (indx = 0; indx < nA; indx++ )
     proclist[ mapA[ indx ] ] = 1;
    
  for (indx = 0; indx < nB; indx++ )
     proclist[ mapB[ indx ] ] = 1;
    
  for (indx = 0; indx < nC; indx++ )
     proclist[ mapC[ indx ] ] = 1;
    
  mproc = 0;
  for (indx = 0; indx < naproc; indx++)
     if( proclist[ indx ] > 0){
       proclist[ mproc ] = indx;
       mproc++;
      }

  for (indx = mproc; indx < naproc; indx++)
     proclist[ indx ] = -1;

  *nproc = mproc;

  return;
}
