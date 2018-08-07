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

void sonenrm_( n, matrixA, mapA, norm, iwork, work, info)
     Integer *n, *mapA, *iwork, *info;
     DoublePrecision *matrixA, *norm, *work;
     
     /*
       
       this subroutine computes the 1-norm
       
       on input:
       
       n = (input) dimension of the matrix A, B
       mapA, mapB, mapZ = (input) integer array of length n holding the processor
       storing the matrices.
       
       matrixA, matrixB, matrixZ = (input) pointers to the appropriate locations
       eval             = (input, unchanged ) array of m DoublePrecision precision eigenvalues
       
       m = number of columns in vecZ
       
       work = the scratch array holding the matrix W;
       requires nvecsW*n + 2*n
       
       output:
       
       W is overwritten by QB.
       
       */
{
  
  Integer ll, k, i, *mapvecA;
  Integer nvecsA;
  Integer me;
  
  DoublePrecision **buff_ptr, **ptrA;
  
  extern Integer mxmynd_();
  
  extern Integer fil_mapvec_();
  extern Integer count_list();
  extern void sonenrm();
  
  /*
    usual story about error handling
    */


  
  ll = *n;
  
  /*
    should perform global sync to check for errors
    */
  
  me = mxmynd_();
  
  buff_ptr = (DoublePrecision **) work ;
  
  mapvecA = iwork;
  nvecsA = fil_mapvec_( &me, &ll, mapA, mapvecA);
  
  ptrA = buff_ptr;
  buff_ptr += nvecsA;
  
  k = 0;
  for ( i = 0; i < nvecsA; i++ )
    {
      ptrA[i] = matrixA + k;
      k += (ll - mapvecA[i]);
    }
  
  sonenrm( n, ptrA, mapA, norm, iwork, work + nvecsA , info);
  
  return;
}

