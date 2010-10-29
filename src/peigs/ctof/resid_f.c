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

void resid_( n, matrixA, mapA, m, matrixZ, mapZ, eval, iwork, work, res, info)
     Integer *n, *mapA, *m, *mapZ, *iwork, *info;
     DoublePrecision *matrixA, *matrixZ, *eval, *work, *res;
     
     /*
       
       FORTRAN wrapper for resid().

       this subroutine computes the residual
       
       res = max_{i} | A z_{i} - \lambda_{i} z_{i} |/( | A | * ulp )
       
       where (lambda, z) is a standard eigen-pair.
       ULP = (machine precision) * (machine base)
       
       |A z_{i} ... | is the infinity-norm,
       |A| is the 1-norm of A,
       res is reasonable if it is of order about 50 or less.
       
       Assume that A is an n x n symmetric matrix
       in packed storage format column-distributed
       across all the participating processors
       
       on input:
       
       n = (input) dimension of the matrix A
       mapA, mapZ = (input) integer array of length n holding the processor
       storing the matrices.
       
       matrixA, matrixZ = (input) pointers to the appropriate locations
       eval             = (input, unchanged ) array of m DoublePrecision precision eigenvalues
       
       m = number of columns in vecZ
       
       work = scratch array 
       */
{
  
  Integer ll, k, i, *iscrat, *mapvecA, *mapvecZ;
  Integer nvecsA, nvecsZ, me;
  
  DoublePrecision *scratch;
  DoublePrecision **buff_ptr, **ptrA, **ptrZ;
  
  
  extern Integer  mxmynd_();
  extern Integer  fil_mapvec_();
  extern void     resid();
  
  /*
    usual story about error handling
    */

  ll = *n;
  
  if ( ll < 1 )
    return;  /* error */
  
  /*
    should perform global sync to check for errors
    */
  
  me = mxmynd_();

  scratch  = work;
  buff_ptr = (DoublePrecision ** ) work;
  
  iscrat = iwork;

  mapvecA = iscrat;
  nvecsA = fil_mapvec_( &me, &ll, mapA, mapvecA);
  iscrat += nvecsA;
  
  mapvecZ = iscrat;
  nvecsZ = fil_mapvec_( &me, m, mapZ, mapvecZ );
  iscrat += nvecsZ;
  
  ptrA = buff_ptr;
  buff_ptr += nvecsA;
  scratch  += nvecsA + 1;
  
  k = 0;
  for ( i = 0; i < nvecsA; i++ ) {
    ptrA[i] = &matrixA[k];
    k += ( ll - mapvecA[i] );
  }
  
  ptrZ = buff_ptr;
  buff_ptr += nvecsZ;
  scratch  += nvecsZ + 1;

  k = 0;
  for ( i = 0; i < nvecsZ; i++ ) {
    ptrZ[i] = &matrixZ[k];
    k += ll;
  }
  
  scratch += 2 * ( nvecsZ + 1 );

  resid( n, ptrA, mapA, m, ptrZ, mapZ,
	 eval, buff_ptr, iwork, scratch, res, info);
  
  return;
}

