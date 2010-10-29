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
#include <stdlib.h>

#include "globalp.c.h"

void ortho_( n, m, matZ, mapZ, iwork, work, ort, info)
     Integer *n, *m, *mapZ, *iwork, *info;
     DoublePrecision *matZ, *work,  *ort;
     
     /*
       
       FORTRAN wrapper for ortho.

       This subroutine computes the infinity-norm measure:
       
       ort = max_i | (Z^t.Z)_i - I_i | / ULP,

       for the standard symmetric eigensystem problem where

         Z is N-by-M
         I is M-by-M.
         Z_i is an eigenvector,
         (Z^t.Z)_i is the i-th column of Z^t.Z
         I_i is the i-th column of the m-by-m identity matrix
         ULP = (machine precision) * (machine base)
         |.| is the infinity-norm.

       res is reasonable if it is of order about 50 or less.


       MUST have M <= N.  If M > N then program exits.
       This is not a limitation of this subroutine as M > N
       implies the columns of Z are linearly dependent, which
       implies "ort" will always be large in this case.

       
       on input:
       
       n    = (input) Integer
              1st dimension of the matrix Z

       m    = (input) Integer
              2nd dimension of the matrix Z

       matZ =(input) double precision
             colZ[i] points to the start of the ith column of Z OWNED
             by this processor.
       
       mapZ = (input) integer array of length m
              mapZ[i] is  the id of the processor owning column i of Z.
       
       iwork = (Workspace) Integer array of length 
       work  = (Workspace) DoublePrecision array of length 
       
       OUTPUT:

       ort   = (Output) DoublePrecision
               ort as defined in equation above, but only if
               info = 0.  
       
       info  = (Output) Integer
               = 0  if this processor is in mapZ
               = -1 if this processor is NOT in mapZ, in which case
                    ort is undefined on exit.
       
       */
{
  
  Integer          i, nvecs, me;

  DoublePrecision *scratch;
  DoublePrecision **colZ, **ibuffptr;
  
  extern Integer mxmynd_();
  extern void    ortho();


  me = mxmynd_();

  if( *n < *m ) {
      fprintf(stderr,
              "Error in routine ortho_. m (=%d) < n (=%d). me = %d \n",
              *m, *n, me);
      exit(-1);
  }
  
  colZ = (DoublePrecision ** ) work;

  nvecs = -1;
  for( i = 0; i < *m; i++ ) {
    if( mapZ[i] == me ) {
       nvecs++;
       colZ[nvecs] = &matZ[*n * nvecs];
    }
  }
  nvecs++;
  
  ibuffptr = (DoublePrecision ** ) work + nvecs + 1;
  scratch  = (DoublePrecision * ) work + 2 * nvecs + 2;


  ortho( n, m, colZ, mapZ, ibuffptr, iwork, scratch, ort, info);
  
  
  return;
}

