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

void bortho_( n, matB, mapB, m, matZ, mapZ, iwork, work, ort, info)
     Integer *n, *mapB, *m, *mapZ, *iwork, *info;
     DoublePrecision *matB, *matZ, *work, *ort;
     
     /*
       
       FORTRAN wrapper for b_ortho.

       This subroutine computes the  infinity-norm measure:
       
       ort = max_i | (Z^t.B.Z)_i - I_i | / ULP,

       for the generalized symmetric eigensystem problem where

         (lambda, z_i ) is a generalized eigen-pair,
         (Z^t.B.Z)_i is the i-th column of Z^t.B.Z
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
              Dimension of the matrix B, and 1st dimension of the matrix Z.

       matB =(input) double precision
             colB[i] points to the start of the ith column of B OWNED
             by this processor.  B is n-by-n symmetric stored in packed
             format.
       
       mapB = (input) integer array of length n
              mapB[i] is  the id of the processor owning column i of B.
       
       m    = (input) Integer
              2nd dimension of the matrix Z

       mapZ = (input) integer array of length m
              mapZ[i] is  the id of the processor owning column i of Z.
       
       matZ =(input) double precision
             colZ[i] points to the start of the ith column of Z OWNED
             by this processor.
       
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
  
  Integer          i, k, nvecsB, nvecsZ, me;

  DoublePrecision *scratch;
  DoublePrecision **colZ, **colB, **ibuffptr;
  
  extern Integer mxmynd_();
  extern void    b_ortho();


  me = mxmynd_();

  if( *n < *m ) {
      fprintf(stderr,
              "Error in routine bortho_. m (=%d) < n (=%d). me = %d \n",
              *m, *n, me);
      exit(-1);
  }
  
/*
 *  Matrix B is triangular
 */

  colB = (DoublePrecision ** ) work;

  k     = 0;
  nvecsB = 0;
  for( i = 0; i < *n; i++ ) {
    if( mapB[i] == me ) {
       nvecsB++;
       colB[nvecsB-1] = &matB[k];
       k += *n - i;
    }
  }
  

/*
 *  Matrix Z is full
 */

  colZ = (DoublePrecision ** ) work + nvecsB + 1;

  nvecsZ = -1;
  for( i = 0; i < *m; i++ ) {
    if( mapZ[i] == me ) {
       nvecsZ++;
       colZ[nvecsZ] = &matZ[*n * nvecsZ];
    }
  }
  nvecsZ++;
  

  ibuffptr = (DoublePrecision ** ) work + nvecsB + nvecsZ + 2;
  scratch  = (DoublePrecision * ) work + nvecsB + 4 * nvecsZ + 5;


  b_ortho( n, colB, mapB, m, colZ, mapZ, ibuffptr, iwork, scratch, ort, info);
  
  
  return;
}

