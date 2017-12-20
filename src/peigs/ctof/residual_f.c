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

#define MSG_START 25000
#define ZERO ((DoublePrecision) 0.0e0)



void residual_( n, matrixA, mapA, matrixB, mapB, m, matrixZ, mapZ, eval, iwork, work, res, info)
     Integer *n, *mapA, *mapB, *m, *mapZ, *iwork, *info;
     DoublePrecision *matrixA, *matrixB, *matrixZ, *eval, *work, *res;
     
     /*
       
       this subroutine computes the residual
       
       res = max_{i} | A z_{i} - \lambda_{i} B z_{i} |/( | A | * ulp )
       
       where (lambda, z) is a generalized eigen-pair.
       ULP = (machine precision) * (machine base)
       
       |A z_{i} ... | is the infinity-norm,
       |A| is the 1-norm of A,
       res is reasonable if it is of order about 50 or less.
       
       
       Assume that A and B are n x n symmetric matrices
       in packed storage format column-distributed
       across all the participating processors
       
       on input:
       
       n = (input) dimension of the matrix A, B
       mapZ = (input) integer array of length m holding the processor
                    storing the matrices.
       mapA, mapB = (input) integer array of length n holding the processor
                    storing the matrices.
       
       matrixA, matrixB, matrixZ = (input) pointers to the appropriate locations
       eval             = (input, unchanged ) array of m DoublePrecision precision eigenvalues
       
       m = number of columns in vecZ
       
       work = scratch array 
       
       */
{
  
  Integer ll, k, i, *iscrat, *mapvecA, *mapvecB, *mapvecZ;
  Integer nvecsA, nvecsB, nvecsZ;
  Integer me;
  
  DoublePrecision *scratch;
  DoublePrecision **buff_ptr, **ptrA, **ptrB, **ptrZ;
  
  
  /*
    blas call
    */
  
  extern Integer mxwrit_(), mxread_();
  extern Integer menode_(), mxmynd_();
  
  extern Integer fil_mapvec_();
  extern Integer count_list();
  extern void residual();
  
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
  
  mapvecA = iwork;
  iscrat = iwork;
  nvecsA = fil_mapvec_( &me, &ll, mapA, mapvecA);
  iscrat += nvecsA;
  
  mapvecB = iscrat;
  nvecsB = fil_mapvec_( &me, &ll, mapB, mapvecB);
  iscrat += nvecsB;
  
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
  
  ptrB = buff_ptr;
  buff_ptr += nvecsB;
  scratch  += nvecsB + 1;
  
  k = 0;
  for ( i = 0; i < nvecsB; i++ ) {
    ptrB[i] = &matrixB[k];
    k += ( ll - mapvecB[i] );
  }
  
  ptrZ = buff_ptr;
  buff_ptr += nvecsZ;
  scratch  += nvecsZ + 1;

  k = 0;
  for ( i = 0; i < nvecsZ; i++ ) {
    ptrZ[i] = &matrixZ[k];
    k += ll;
  }
  
  scratch  += 3 * ( nvecsZ + 1 );

  residual( n, ptrA, mapA, ptrB, mapB, m, ptrZ, mapZ,
	    eval, buff_ptr, iwork, scratch, res, info);
  
  return;
}

