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

void tresid_( n, m, d, e, matrixZ, mapZ, eval, iwork, work, res, info)
     Integer *n, *m, *mapZ, *iwork, *info;
     DoublePrecision *d, *e, *matrixZ, *eval, *work, *res;
     
     /*
       
       FORTRAN wrapper for tresid().

       this subroutine computes the residual
       
       res = max_{i} | T z_{i} - \lambda_{i} z_{i} |/( | T | * ulp )
       
       where T is an n-by-n  tridiagonal matrix,
        ( \lambda_{i} , z_{i} ) is a standard eigen-pair, and
       ULP = (machine precision) * (machine base)
       
       |T| is the 1-norm of T
       |T z_{i} .... | is the infinity-norm
       res is reasonable if it is of order about 50 or less.
       
   Arguments
   ---------
       
       n = (input) dimension of the matrix T
       m = (input) number of eigenvalues/eigenvectors

       d(1:n) = (input) diagonal of T

       e(1)   = junk
       e(2:n) = (input) sub-diagonal of T
                        = super-diagonal of T

       colZ = (input) DoublePrecision eigenvectors in packed storage
       mapZ = (input) integer array of length m holding the id of the processors
                      holding the eigenvectors in colZ.
       eval = (input ) array of m DoublePrecision eigenvalues
       
       iwork = (workspace) array of Integer workspace
       work = (workspace) array of DoublePrecision workspace
       
       res  = (output ) the residual described above.

       info = (output ) = 0, then everything ok
                        < 0, then -info-th input variable is not valid


       */
{
  
  Integer ll, k, i;
  Integer nvecsZ, me;
  
  DoublePrecision *scratch;
  DoublePrecision **buff_ptr, **ptrZ;
  
  
  extern Integer  mxmynd_();
  extern Integer  count_list();
  extern void     resid();
  
  /*
    usual story about error handling
    */

  ll = *n;
  
  if ( ll < 1 )
    return;  /* error */
  
  me = mxmynd_();

  scratch  = work;
  buff_ptr = (DoublePrecision ** ) work;
  
  nvecsZ = count_list( me, mapZ, m );
  
  ptrZ = buff_ptr;
  buff_ptr += nvecsZ + 1;
  scratch  += nvecsZ + 1;

  k = 0;
  for ( i = 0; i < nvecsZ; i++ ) {
    ptrZ[i] = &matrixZ[k];
    k += ll;
  }
  
  tresid( n, m, d, e, ptrZ, mapZ, eval, iwork, scratch, res, info);
  
  return;
}

