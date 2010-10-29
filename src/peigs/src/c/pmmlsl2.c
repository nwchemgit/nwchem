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
#include <math.h>

#include "globalp.c.h"

void pmmlsl2(n, mapA, mapvecA, vecA, mapB, mapvecB, vecB,
	     iscratch, dscratch, buff_ptr)
     
     Integer            *n, *mapA, *mapvecA, *mapB, *mapvecB, *iscratch;
     DoublePrecision         *dscratch;
     DoublePrecision        **vecA, **vecB, **buff_ptr;
{
  /*
   * ------------------------------------------------------------------------
   *
   * This subroutine computes the lower triangular part of the product of a
   * lower triangular matrix, B, with a symmetric matrix, A, and
   * stores the result in A.  It is assumed that both A and B are square
   * matrices.
   *
   * A <- lower triangular part of B * A
   *
   * NOTE: B * A is not a symmetric matrix, thus computing only the lower
   *       triangular portion of B * A, as this routine does, is not the
   *       same as computing all of B * A.
   *
   *
   * A: a symmetric matrix with the lower triangular portions of its columns
   *    stored in distributed column format (or equivalently the upper
   *    triangular portions of its rows stored in distributed column format).
   *
   * B: a lower triangular matrix stored in distributed column format.
   *
   * n         = dimension of the matrix
   * mapA      = index array holding the proces
   * mapvecA   = index array holding the i-th column this processor owns
   * vecA      = DoublePrecision pointer to the array location of the i-th column of A
   * mapB      = index array holding the proces
   * mapvecB   = index array holding the i-th column this processor owns
   * vecB      = DoublePrecision pointer to the array location of the i-th column of B
   *
   * dscratch  = (workspace) (DoublePrecision *)
   * iscratch  = (workspace) Integer size 
   * buff_ptr  = (workspace) array of pointers to DoublePrecision precision numbers;
   *                         size 
   *
   * ------------------------------------------------------------------------
   */
  
  /*
   * Local Variables
   * ---------------
   */

  static Integer      ONE = 1;

  Integer             me, isize, indx, nvecsA, nvecsB;
  DoublePrecision *d_ptr, *scratch, *buff, **b_ptr, **matrix2;

  /*
   * External procedures
   * -------------------
   */

  extern void     dcopy_();
  extern Integer      mxmynd_(), mxwrit_(), mxread_();
  extern void     zero_out(), de_sym(), mxm_ll();

  /*
   * Executable code
   * ---------------
   */

  me = mxmynd_();

  /*
   *  Partition workspace.
   */
  
  scratch = dscratch;
  
  nvecsA = count_list(me, mapA, n);
  nvecsB = count_list(me, mapB, n);
  
  /*
   * Copy the A matrix over to "buff"
   */
  
#ifdef DEBUG
  fprintf(stderr, " just in pmmlsl2 1 \n");
#endif
  
  b_ptr = matrix2 = buff_ptr;
  b_ptr +=  nvecsA;
  buff    = dscratch;
  
  /*
    first step in un-symmetrizing matrix A
    copy the lower triangular part into a full matrix
    */

  for ( indx = 0; indx < nvecsA; indx++ ){
    matrix2[indx] = buff + *n*indx;
    d_ptr = matrix2[indx];
    zero_out(*n, d_ptr);
    isize  = *n - mapvecA[indx];
    d_ptr += mapvecA[indx];
    dcopy_(&isize, vecA[indx], &ONE, d_ptr, &ONE);
  }
  buff += nvecsA * *n;
  
#ifdef DEBUG
  fprintf(stderr, " just in pmmlsl2 2 \n");
#endif
  
  /*
    all to all swap for explicitely representing the symmetric matrix
    */
  
#ifdef DEBUG1
  fprintf(stderr, " just before de_sym, me = %d  \n", me);
#endif
  
  de_sym( *n, 1111, buff, vecA, mapA, matrix2, iscratch, b_ptr);

#ifdef DEBUG1
  fprintf(stderr, " after de_sym, me = %d  \n", me);
#endif
  
  
  /*
    matrix2 is the full symmetric matrix of vecA
    multiply_lower_triangular times full matrix
    */
  
#ifdef DEBUG  
  fprintf(stderr, " ***************************** out of de_sym just before mxm_llx, me = %d  \n", me);
#endif

#ifdef DEBUG1
  fprintf(stderr, " before mxm_llx, me = %d ***************  \n", me);
#endif
  
  mxm_llx( n, vecB, mapB, n, matrix2, mapA, iscratch, buff);
  
#ifdef DEBUG1
  fprintf(stderr, " after mxm_llx, me = %d ******************   \n", me);
#endif
  
  
  /*
    copy back
    */
  
  for ( indx = 0; indx < nvecsA; indx++ ){
    isize = *n - mapvecA[indx];
    dcopy_(&isize, matrix2[indx] + mapvecA[indx], &ONE, vecA[indx], &ONE);
  }
  
  return;
}

