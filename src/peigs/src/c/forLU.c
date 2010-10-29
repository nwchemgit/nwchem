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

/*
  
  PeIGS internal routine.  This is used in the one processor case.
  The user should be warned of possible communication hang in the
  multiprocessor case if there's not enough buffer space.
  
  
  Let U be any upper triangular matrix. Let L be a lower triangular matrix.
  This subroutine solves the upper-triangular part of Y for the matrix
  problem
  
  L Y = U
  
  The resulting upper-triangular part of matrix Y overwrites U.
  
  The matrix U is row wrapped.
  The resulting matrix Y overwrites the matrix U.
  The matrix L is column wrapped.
  
  */

void forwardLU_ ( n, mapL, mapvecL, colL, mapU, mapvecU, rowU, buffer, nprocs, proclist, ibuffer, buff, buff_ptr)
     Integer *n, *mapL, *mapvecL, *mapU, *mapvecU, *nprocs, *proclist, *ibuffer;
     DoublePrecision **colL, **rowU, *buffer, *buff, **buff_ptr;
{
  /*
    
    n            = dimension of the matrix
    mapL         = index array holding the proces
    mapvecL      = index array holding the i-th column this processor owns
    colL         = DoublePrecision pointer to the array location of the i-th column
    
    ditto for mapvecU, rowU, mapU
    buffer      = (DoublePrecision ) of size 2n, buffer for message passing
    ibuffer     = integer buffer space of size p
    
    remark: one can improve the x-fer by taking advantage size of the vector
    
    */
  
  static Integer ONE = 1, MINUSONE = 99999;
  Integer i, k, me, isize, indxx;
  Integer i_L, i_U, nvecsU;
  
  DoublePrecision t, *d_ptr;
  extern Integer count_list();
  extern void pipe_bcst_fut_col();
  
  /*
    mxsubs calls
    */
  
  extern Integer mxmynd_ ();
  
  /*
    blas calls
    */
  
  extern void dcopy_ ();          
  
  /* 
    mxsubs calls
    */
  
  extern Integer mxwrit_ (),  mxread_ ();
  extern void dscal_();
  extern void daxpy_();
  
  me = mxmynd_ ();
  
  nvecsU = count_list ( me, mapU, n);
  i_L=0;
  i_U=0;
  
  for ( i = 0; i < *n; i++ ) {   /* indexing on L[i] */
    if ( mapL[i]== me ) {
      if ( mapU[i] != me ) {
	isize = *n - i;
	isize = isize*sizeof(DoublePrecision);
	isize = mxwrit_ ( colL[i_L], &isize, &mapU[i], &MINUSONE);
	isize = 2 * ( *n -i ) *sizeof(DoublePrecision);
	pipe_bcst_fut_col( *n , (char *) buffer, isize, i, mapU[i], i, mapU, ibuffer);
      }
      
      if ( mapU[i] == me ){
	/*
	  i own both U[i] and L[i]; copy them into one vector
	  buffer buffer owns L[i:n-1] U[n-1:i]
	  this is a trick to minimize data movement but only pass the 
	  necessary vectors
	  */
	
	t = (DoublePrecision) 1.0e0 / (DoublePrecision) *colL[i_L] ;
	isize = *n - i;
	d_ptr = rowU[i_U];
	dscal_(&isize, &t, d_ptr, &ONE);
	d_ptr = colL[i_L];
	dcopy_( &isize,  d_ptr, &ONE, buffer, &ONE );
	d_ptr = rowU[i_U];
	dcopy_( &isize,  d_ptr, &ONE, &buffer[isize] , &ONE );
	
	/*
	  broadcasting l[i,i]... l[n-1, i] y[i,i] ... y[n-1 , i ] ; just a choice since
	  a better way is to map it via l[i,i], ... l[n-1, i], y[n-1, i], ... y[i,i]
	  it removes the problems of shifting later
	  */
	
	isize = 2*(*n - i)*sizeof(DoublePrecision);
	pipe_bcst_fut_col( *n , (char *) buffer, isize, i, mapU[i], i, mapU, ibuffer);
	i_U++;
      }
      i_L++;	  
    }
    
    if ( mapL[i] != me ) {
      if ( mapU[i] == me ){
	isize = (*n - i)*sizeof(DoublePrecision);
	isize = mxread_ ( buffer, &isize, &mapL[i], &MINUSONE);
	t = *buffer; /* this is L[i,i] */
	t = (DoublePrecision) 1.0e0 / t;
	isize = *n - i;
	d_ptr = rowU[i_U];
	dscal_(&isize, &t, d_ptr, &ONE);
	dcopy_( &isize,  d_ptr, &ONE, &buffer[isize], &ONE );
	isize = 2*(*n - i)*sizeof(DoublePrecision);
	pipe_bcst_fut_col( *n , (char *) buffer, isize, i, mapU[i], i, mapU, ibuffer);
	i_U++;
      }
      else
	if ( mapU[i] != me ) {
	  /* receive the Ls and the Us */
	  /*
	    if ( i != 0 ){
	    */
	  isize = 2*(*n-i)*sizeof(DoublePrecision); 
	  pipe_bcst_fut_col( *n , (char *) buffer, isize, i, mapU[i], i, mapU, ibuffer );
	}
    }
    
    /*
      buffer now contains the ( L[i:n-1], Y[i:n-1] ); this is the vectors that is coming in
      */
    
    for ( k = i_U; k < nvecsU; k++ ){
      indxx = mapvecU[k];
      t = (DoublePrecision) -1.0e0 * buffer[indxx - i];
      
      /*
	this is u(i, indx)
	*/
      
      isize = *n - indxx ;
      d_ptr = rowU[k];
      daxpy_ ( &isize, &t, &buffer[*n - i + indxx - i] , &ONE, d_ptr, &ONE);
    }
  }
  return;
}
