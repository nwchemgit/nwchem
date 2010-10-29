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

void pmmLSL(n, mapA, mapvecA, vecA, mapB, mapvecB, vecB,
            iscratch, dscratch, buff_ptr)
     
     Integer            *n, mapA[], mapvecA[], mapB[], mapvecB[], iscratch[];
     
     DoublePrecision         dscratch[];
     
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
   * dscratch  = (workspace) (DoublePrecision *) size (2n + size of A)
   * iscratch  = (workspace) Integer size (p = # of distinct processor ids in mapA)
   * buff_ptr  = (workspace) array of pointers to DoublePrecision precision numbers;
   *                         size (same size as vecA = nvecsA)
   *
   * ------------------------------------------------------------------------
   */
  
  /*
   * Local Variables
   * ---------------
   */
  
  static Integer      ONE = 1, MINUSONE = -1;
  
  Integer             i, k, me, isize, indx, i_A, i_B, nvecsA, nvecsB;
  Integer             num_procs;
  
  DoublePrecision         *d_ptr, *t_ptr, *scratch, *buff;
  
  
  /*
   * External procedures
   * -------------------
   */
  
  extern void     dcopy_();
  extern Integer      mxmynd_(), mxwrit_(), mxread_();
  extern Integer      count_list();
  extern void         chol_pipe_bcast();
  extern Integer reduce_list2();
  extern void daxpy_();
  extern void     zero_out();
  
  /*
   * Executable code
   * ---------------
   */
  
  me = mxmynd_();
  
  /*
   *  Partition workspace.
   */
  
  scratch = dscratch;
  buff    = dscratch + 2 * *n;
  
  nvecsA = count_list(me, mapA, n);
  nvecsB = count_list(me, mapB, n);
    
  /*
   * Copy the A matrix over to "buff"
   */
  
  d_ptr = buff;
  for (i = 0; i < nvecsA; i++){
    isize       = *n - mapvecA[i];
    buff_ptr[i] = d_ptr;
    dcopy_(&isize, vecA[i], &ONE, d_ptr, &ONE);
    d_ptr += isize;
    zero_out(isize, vecA[i]);
  }
  
  i_A = 0;
  i_B = 0;
  num_procs = reduce_list2( *n, mapA, iscratch );
  
  for (i = 0; i < *n; i++){
    if (mapA[i] == me){
      /*
       * Copy A[i:n-1] B[n-1:i] to scratch.
       */
      
      if (mapB[i] == me){
	isize = *n - i;
	dcopy_(&isize, buff_ptr[i_A], &ONE, scratch,         &ONE);
	dcopy_(&isize, vecB[i_B],     &ONE, scratch + isize, &MINUSONE);
	isize *= 2 * sizeof(DoublePrecision);
	/*
	  chol_pipe_bcast(scratch, isize, i, mapA[i], *n, mapA, iscratch);
	  */
	bbcast00( (char *) scratch, isize, i, mapA[i], num_procs, iscratch );
	i_B++;
      }
      else {
	isize = (*n - i) * sizeof(DoublePrecision);         
	isize = mxread_(scratch, &isize, &mapB[i], &isize);
	
	isize = *n - i;
	dcopy_(&isize, scratch,       &ONE, scratch + isize, &MINUSONE);
	dcopy_(&isize, buff_ptr[i_A], &ONE, scratch,         &ONE);
	isize = 2 * (*n - i) * sizeof(DoublePrecision);
	/*
	  chol_pipe_bcast(scratch, isize, i, mapA[i], *n, mapA, iscratch);
	  */
	bbcast00( (char *) scratch, isize, i, mapA[i], num_procs, iscratch );
      }
      i_A++;
    }
    else
      {
	if (mapB[i] == me) {
	  isize = (*n - i) * sizeof(DoublePrecision);
	  isize = mxwrit_(vecB[i_B], &isize, &mapA[i], &isize);
	  i_B++;
	}
	
	isize = 2 * (*n - i) * sizeof(DoublePrecision);
	/*
	  chol_pipe_bcast(scratch, isize, i, mapA[i], *n, mapA, iscratch);
	  */
	bbcast00( (char *) scratch, isize, i, mapA[i], num_procs, iscratch );
      }
    
    /*
     * Scratch now contains A[i:n-1] with B[n-1:i]
     */
    
    d_ptr = scratch + *n - i;         /* points to  b[i, i] */
    for (k = 0; k < nvecsA; k++) {
      indx = mapvecA[k];
      if (indx >= i){
	isize = *n - indx;
	t_ptr = scratch + indx - i;     /* points to a[indx, i+(indx-i)] */
	daxpy_(&isize, t_ptr, d_ptr, &MINUSONE, vecA[k], &ONE);
      }
      else{
	isize = *n - i;
	t_ptr = buff_ptr[k] + i - indx;  /* points to a[i+(i-indx), indx] */
	daxpy_(&isize, t_ptr, d_ptr, &MINUSONE, vecA[k] + i - indx, &ONE);
      }
    }
  }
  return;
}

void zero_out(n, array)
     Integer             n;
     DoublePrecision         *array;
{
  
  static DoublePrecision      ZERO = 0.0e0;
  Integer             i;
  DoublePrecision         *ptr;
  
  ptr = array;
  for (i = 0; i < n; i++)
    *ptr++ = ZERO;
  
  return;
}
