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
/*
  PeIGS internal routine forwardLL
  
  Let L, Z be two non-singular lower triangular matrices. Let Y be a
  symmetric matrix.  This subroutine solves the lower-triangular part of
  Y for the matrix problem
  
  L Y = Z
  
  The resulting lower triangular part of matrix Y overwrites Z.
  
  Data structure:
  
  The matrix L and Z are assumed to be column wrapped.
  
  !!!! The user should be warned that this routine will hang if
  insufficient message buffering is not available.  PeIGS uses this
  routine for the one-processor case. !!!!!

  */

#include <stdio.h>
#include <math.h>

#include "globalp.c.h"


void forwardLL_ ( n, mapL, mapv_L, colL, mapZ, mapvecZ, colZ, scratch, nprocs, proclist, iscratch)
     Integer *n, *mapL, *mapv_L, *mapZ, *mapvecZ, *nprocs, *proclist, *iscratch;
     DoublePrecision **colL, **colZ, *scratch;
{
  /*
    
    n            = dimension of the matrix
    mapL         = index array holding the proces
    mapv_L      = index array holding the i-th column this processor owns
    colL         = DoublePrecision pointer to the array location of the i-th column
    
    ditto for mapvecU, rowU, mapZ
    scratch      = (DoublePrecision ) of size 2n, scratch buffer for message passing
    iscratch     = integer scratch space of size p
    
    remark: one can improve the x-fer by taking advantage size of the vector by segmenting the
            column vectors

	    this is a communication intensive version since I can't figure out how
	    to code a more elegant version ( 4/1/93 gif)
	    
	    
*/

  static Integer ONE = 1;
  Integer me, isize, indx, indx1, itype;
  Integer iL, iZ, nvecsL, nvecsZ;
  Integer ii, indx2;
  DoublePrecision t;
  
  extern Integer mxmynd_ ();
  extern void combine_vector();
  
  /* blas calls */
  
  extern void dcopy_ ();          
  
  /* 
    mxsubs calls
    */
  
  extern Integer mxwrit_ (),  mxread_ (), indxL();
  extern Integer count_list();
  extern void daxpy_();
  
  me = mxmynd_ ();

  
  nvecsZ = count_list(me, mapZ, n);
  nvecsL = count_list(me, mapL, n);

  if (nvecsZ == *n  && nvecsL == *n && mapZ[0] == mapL[0])
  {
    *nprocs = 1;
  }
  else
  {
    *nprocs = 2;
  }

/*
  fprintf(stderr, " nvecsZ = %d nvecsL = %d \n", nvecsZ, nvecsL);
*/
  
  for ( iZ = 0; iZ < *n; iZ++ ) {
    for ( iL = 0; iL < *n; iL++ ) {
      /*
      fprintf(stderr, " me = %d  iZ = %d  iL = %d \n", me, iZ, iL );
	first case
	*/
      if ( iZ == iL ) {
	if ( mapZ[iZ] != me ) {
	  for ( ii = 0; ii < iZ+1; ii++ ) {
	    if ( mapL[ii] == me ) { 
	      itype = *n + ii;
	      indx = indxL ( ii, nvecsL, mapv_L );
	      isize = (*n - iZ)*sizeof(DoublePrecision);
	      isize = mxwrit_ ( colL[indx] + iZ - ii , &isize, &mapZ[iZ] , &itype);
	    }
	  }
	}
	
	
	if ( mapZ[iZ] == me ) {
	  indx2 = indxL ( iZ, nvecsZ, mapvecZ );
	  if ( *nprocs  > 1 ) {
	    for ( ii = 0; ii < iZ; ii++ ) {
	      if ( mapL[ii] == me ) {
		/*
		  i own both the L vector and the Z vector
		  */
		indx = indxL ( ii, nvecsL, mapv_L );
		isize = *n - iZ;
		dcopy_ ( &isize, colL[indx] + iZ - ii , &ONE, scratch, &ONE);
	      }
	      else {
		isize = (*n - iZ) * sizeof(DoublePrecision);
		itype = *n + ii;
		isize = mxread_ ( scratch, &isize, &mapL[ii], &itype);
	      }
	      
	      
	      if ( mapZ[ii] == me ) {
		indx1 = indxL ( ii, nvecsZ, mapvecZ );
		t = - colZ[indx1][iZ-ii];
	      }
	      else {
		itype = iL;
		isize = sizeof(DoublePrecision);
		isize = mxread_ ( &t , &isize, &mapZ[ii] , &itype);
	      }
	      /*
		scratch contains the L vectors ; L[iL : n-1, iL]
		*/
	      isize = *n - iZ;	      
	      daxpy_ ( &isize, &t, scratch, &ONE, colZ[indx2], &ONE);
	    }
	    
	    /*
	      the ii = iZ case
	      */
	    if ( mapL[iL] == me ) {
	      /*
		i own both the L vector and the Z vector
		*/
	      indx = indxL ( iL, nvecsL, mapv_L );
	      isize = *n - iL;
	      dcopy_ ( &isize, colL[indx], &ONE, scratch, &ONE);
	    }
	    else {
	      isize = (*n - iZ) * sizeof(DoublePrecision);
	      itype = *n + iZ;
	      isize = mxread_ ( scratch, &isize, &mapL[iL], &itype);
	    }
	    
	    indx = indxL ( iZ, nvecsZ, mapvecZ );
	    *colZ[indx] = *colZ[indx]/ *scratch;
	    t = - *colZ[indx];
	    isize = *n - iZ - 1;
	    daxpy_ ( &isize, &t, scratch + 1, &ONE, colZ[indx2] + 1, &ONE);
	  }
	  
	  if ( *nprocs == 1 ) {
	    indx = indxL ( iZ, nvecsZ, mapvecZ );
	    isize = *n - iZ;
	    for ( ii = 0; ii < iZ; ii++ ) {
	      indx1 = indxL ( ii, nvecsL, &mapv_L[0] );
	      dcopy_ ( &isize, colL[indx1] + iZ - ii, &ONE, scratch, &ONE);
	      t = - colZ[ii][iZ-ii];
	      daxpy_ ( &isize, &t, scratch, &ONE, colZ[indx], &ONE);
	    }
	    
	    /*
	      the case iL = iZ
	      */
	    
	    indx1 = indxL ( iL, nvecsL, mapv_L );
	    dcopy_ ( &isize, colL[indx1], &ONE, scratch, &ONE);
	    indx = indxL ( iZ, nvecsZ, mapvecZ );
	    *colZ[indx] = *colZ[indx]/ *scratch;
	    t = - *colZ[indx];
	    isize = *n - iZ - 1;
	    daxpy_ ( &isize, &t, scratch + 1, &ONE, colZ[indx] + 1, &ONE);
	  }
	}
      }
      
      if ( iL > iZ ) {
	if ( mapZ[iZ] == me ) {
	  if ( mapL[iL] == me ) {
	    isize = *n - iL;
	    indx = indxL ( iL, nvecsL, mapv_L );
	    dcopy_ ( &isize, colL[indx], &ONE, scratch, &ONE);
	  }
	  else {
	    isize = *n - iL;
	    isize *= sizeof(DoublePrecision);
	    itype = *n + iL;
	    isize = mxread_ ( scratch, &isize, &mapL[iL], &itype);
	  }
	  
	  /*
	    scratch contains the L[iL][iL:n-1] vector
	    */
	  
	  indx = indxL ( iZ, nvecsZ, mapvecZ ); /* storage location of this column */
	  itype = iL - iZ;
	  colZ[indx][itype] = colZ[indx][itype] / *scratch;
	  t = - colZ[indx][itype];
	  
	  if ( mapZ[iL] != me ){
	    itype = iL;
	    isize = sizeof(DoublePrecision);
	    isize = mxwrit_ ( &t, &isize, &mapZ[iL], &itype);
	  }
	  
	  isize = *n - iL - 1;
	  daxpy_ ( &isize, &t, scratch + 1, &ONE, colZ[indx] + iL - iZ + 1, &ONE );
	}
	else {
	  if ( mapL[iL] == me ) { 
	    isize = *n - iL;
	    isize  = isize * sizeof( DoublePrecision );
	    indx = indxL ( iL, nvecsL, mapv_L );
	    itype = *n + iL;
	    isize = mxwrit_ ( colL[indx] , &isize, &mapZ[iZ] , &itype);
	  }
	}
      }
    }
  }
  return;
}

