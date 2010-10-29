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

void pmmLUL(n, mapA, mapvecA, vecA, mapB, mapvecB, vecB,
            iscratch, dscratch, buff_ptr)

   Integer            *n, *mapA, *mapvecA, *mapB, *mapvecB, *iscratch;

   DoublePrecision         *dscratch;

   DoublePrecision        **vecA, **vecB, **buff_ptr;

{

/*
 * ------------------------------------------------------------------------
 * This subroutine computes the lower triangular part of the product of a
 * lower triangular matrix, A, with an upper triangular matrix, B, and
 * stores the result in A.  It is assumed that both A and B are square
 * matrices.
 *
 * A <- lower triangular part of A * B
 *
 * NOTE: A * B is not a symmetric matrix, thus computing only the lower
 *       triangular portion of A * B, as this routine does, is not the
 *       same as computing all of A * B.
 *
 *
 * A: a symmetric matrix with the lower triangular portions of its columns
 *    stored in distributed column format (or equivalently the upper
 *    triangular portions of its rows stored in distributed column format).
 *
 * A: a lower triangular matrix with the lower triangular portions of its
 *    columns stored in distributed column format.
 * B: an upper triangular matrix with the upper triangular portions of its
 *    rows stored in distributed rows format.
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

   Integer             i, k, me, isize, isize0, indx, i_A, i_B, nvecsA, nvecsB;
   Integer             num_procs;

   DoublePrecision          t;

   DoublePrecision         *d_ptr, *t_ptr, *scratch, *buff;

/*
 * External procedures
 * -------------------
 */

   extern void     dcopy_();

   extern Integer      mxmynd_(), mxwrit_(), mxread_();
   extern Integer count_list();
   extern Integer reduce_list2();
   extern void chol_pipe_bcast();
   extern void     zero_out();

   extern void daxpy_();

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
   num_procs =reduce_list2( *n, mapA, iscratch );

   nvecsA = count_list(me, mapA, n);
   nvecsB = count_list(me, mapB, n);

   /*
    * Copy the A matrix over to "buff"
    */

   d_ptr = buff;
   for (i = 0; i < nvecsA; i++)
     {
       isize       = *n - mapvecA[i];
       buff_ptr[i] = d_ptr;
       
       dcopy_(&isize, vecA[i], &ONE, d_ptr, &ONE);
       zero_out(isize, vecA[i]);
       
       d_ptr += isize;
     }
   
   i_A = 0;
   i_B = 0;

   for (i = 0; i < *n; i++)
     {
      if (mapA[i] == me)
	{
         /*
          * Copy A[i+1:n-1] B[n-1:i+1] to scratch and process column i.
          */
	  
	  if (mapB[i] == me)
	    {
	      
	      isize = *n - i - 1;
              if( isize > 0 ) {
	        dcopy_(&isize, vecB[i_B] + 1,     &ONE, scratch+isize, &MINUSONE);
	        dcopy_(&isize, buff_ptr[i_A] + 1, &ONE, scratch,         &ONE);
              }
	      
	      isize = 2 * (*n - i - 1) * sizeof(DoublePrecision);
	      /*
		chol_pipe_bcast(scratch, isize, i, mapA[i], *n-i, &mapA[i], iscratch);
		*/
	      
	      bbcast00( (char *) scratch, isize, i, mapA[i], num_procs, iscratch );
	      t     = *vecB[i_B];
	      isize = *n - i;
	      daxpy_(&isize, &t, buff_ptr[i_A], &ONE, vecA[i_A], &ONE);
	      
	      i_B++;
	      
	    }
	  else
	    {
	      
	      isize = (*n - i) * sizeof(DoublePrecision);
	      isize = mxread_(scratch, &isize, &mapB[i], &isize);

	      t = *scratch;              /* = u[i,i] */
	      
	      /*
	       *  B[i][n] is already in the right place in scratch so don't
	       *  move it.  It wastes time (and would also cause dcopy to
	       *  modified an "aliased" variable (which is illegal in FORTRAN).
	       */
	      
	      isize  = *n - i - 1;
	      isize0 = isize-1;
              if( isize0 > 0 )
	         dcopy_(&isize0, scratch + 1,       &ONE, scratch+isize+1, &MINUSONE);
              if( isize > 0 )
	         dcopy_(&isize,  buff_ptr[i_A] + 1, &ONE, scratch, &ONE);
	      
	      isize = 2 * (*n - i - 1) * sizeof(DoublePrecision);
	      /*
		chol_pipe_bcast(scratch, isize, i, mapA[i], *n-i, &mapA[i], iscratch);
		*/
	      
	      bbcast00( (char *) scratch, isize, i, mapA[i], num_procs, iscratch );
	      

	      isize = *n - i;
	      daxpy_(&isize, &t, buff_ptr[i_A], &ONE, vecA[i_A], &ONE);
	      
	    }
	  i_A++;
	}
      else
	{
	  
	  if (mapB[i] == me)
	    {
	      isize = (*n - i) * sizeof(DoublePrecision);
	      isize = mxwrit_(vecB[i_B], &isize, &mapA[i], &isize);
	      i_B++;
	    }
	  
	  isize = 2 * (*n - i - 1) * sizeof(DoublePrecision);
	  /*
	    chol_pipe_bcast(scratch, isize, i, mapA[i], *n-i, &mapA[i], iscratch);
	    */
	  bbcast00( (char *) scratch, isize, i, mapA[i], num_procs, iscratch );
	}
      
      /*
       * Scratch now contains A[i+1:n-1] with B[n-1:i+1]
       */

      d_ptr = scratch + 2 * (*n - i - 1) - 1;                /* = u(i,i+1) */
      for (k = 0; k < nvecsA; k++)
      {
         indx = mapvecA[k];
         if (indx > i)
         {
            t_ptr = d_ptr - (indx - (i + 1));             /* u(i, indx) */
            isize = *n - indx;
            daxpy_(&isize, t_ptr, scratch + indx - i - 1, &ONE, vecA[k], &ONE);
         }
      }
   }

   return;
}




