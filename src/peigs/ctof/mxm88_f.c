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
#define min(a,b) ((a) < (b) ? (a) : (b))

#include "globalp.c.h"

void mxm88_( n, matA, mapA, m, matB, mapB, iwork, work)
     Integer          *n, *mapA, *m, *mapB, *iwork;
     DoublePrecision  *matA, *matB, *work;
{
  /*
   *
   * This is the FORTRAN wrapper for the matrix multiplication
   * routine: mxm88.
   *
   */
  
  Integer me, nvecsA, nvecsB;
  Integer i, k;

  DoublePrecision **iptr;
  
  DoublePrecision **ptrA, **ptrB;
  DoublePrecision *scratch;
  
  extern Integer mxmynd_();
  extern Integer count_list();

  extern void mxm88();



  me = mxmynd_();

  nvecsA = count_list( me, mapA, n);
  nvecsB = count_list( me, mapB, m);
  
  ptrA =(DoublePrecision **) work;
  ptrB =(DoublePrecision **)( work + nvecsA + 1);
  iptr =(DoublePrecision **)( work + nvecsA + nvecsB + 2);
  scratch =(DoublePrecision *)( work + nvecsA + nvecsB + nvecsB + 3 );
  
/*
 *  Matrix A in A * B is triangular or symmetric.
 */

  k     = 0;
  nvecsA = 0;
  for( i = 0; i < *n; i++ ) {
    if( mapA[i] == me ) {
       nvecsA++;
       ptrA[nvecsA-1] = &matA[k];
       k += *n - i;
    }
  }
  
/*
 *  Matrix B in A * B is full.
 */

  for( i = 0; i < nvecsB; i++ ) {
    ptrB[i] = &matB[*n * i];
  }
  

  mxm88( n, ptrA, mapA, m, ptrB, mapB, iwork, scratch, iptr );
  
  return;
}
