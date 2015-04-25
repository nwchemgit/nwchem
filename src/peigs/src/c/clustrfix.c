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
 *  -- PEIGS  routine (version 3.0) --
 *     Pacific Northwest Laboratory
 *     July 28, 1995
 *
 *======================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "globalp.c.h"

#define min(a,b) (((a)) < ((b)) ? (a) : (b))
#define max(a,b) (((a)) > ((b)) ? (a) : (b))
#define ffabs(a) ((a) >= ((DoublePrecision) 0.e0) ? (a) : (-a))


#define R_ZERO (DoublePrecision) 0.0e0
#define R_ONE  (DoublePrecision) 1.0e0
#define R_TEN  (DoublePrecision) 1.0e1
#define ODM3 (DoublePrecision) 1.0e-3
#define ODM1 (DoublePrecision) 1.0e-1
#define EXTRA 1
#define KK    2
#define MAXITS 5
#define DOUBLE 200
#define INTEGER 20
#define TWO  2

#define I_ZERO 0

static Integer clustr_check(Integer, Integer, Integer, Integer);

static Integer clustr_check(c1, cn, imin, imax)
     Integer c1, cn, imin, imax;
{
  /*
     routine to determine if
     the cluster is
     actually in the desired region of the clustr finder
     */
  
  if ( cn < imin )
    return(-1);
  
  if ( imax  < c1 )
    return(-1);
  
  return(1);
}

static DoublePrecision relgap(l1, l2)
     DoublePrecision l1, l2;
{
  
  /*
     assume l1 & l2 != 0
     */
  
  DoublePrecision alpha, beta;
  
  alpha=ffabs(l1-l2);
  beta = alpha/(max(ffabs(l1), ffabs(l2)));
  
  return(beta);
}

Integer clustrfix_ (n, d, e, m, w, iblock, nsplit, isplit, num_clustr, clustr_info )
     Integer n, m, *iblock, nsplit, *isplit, *num_clustr, *clustr_info;
     DoublePrecision *d, *e, *w;
     /*
       this routine finds the cluster information for symmetric tridiagonal matrices
       
       on input:
       
       n = dimension of the matrix
       d = array of n DoublePrecision; diagonal of the tridiagonal matrices
       e = array of n DoublePrecision; e[0] is junk -- left over from eispack--should remove
       super-diagonal of the tridiagonal matrices; 
       m = number of the number eigenvalues; from dstebz
       w = array of m DoublePrecision precision; eigenvalues from dstebz
       iblock = array of n integers ( at most ); from dstebz
       nsplit = number of blocks; output of dstebz
       isplit = split points of blocks; output of dstebz
       ptb_eval = perturbed eigenvalues
       
       output:
       
       function returns the maximum cluster size ;
       
       num_one_blk = number of size 1 blocks
       one_block   = integer array ; location of size 1 blocks
       clustr_info = integer array ; max size 4 * n; cluster information : beg block, end of block
       beg of cluster, end of cluster
       nsplit      = number of clusters

       */
{
  Integer jblk, nblk;
  Integer j, num_cls, num_all_cls;
  Integer b1, num_eig;
  Integer max_clustr_size;
  Integer beg_of_block, end_of_block;
  Integer bn;
  Integer blksiz, c1=-1, c2=-1, ptr;
  Integer me;
  Integer *c_ptr;
  Integer tail;
  Integer clustr_size=-1, iii, jjj;
  
  DoublePrecision eps, xj1, xj2;
  DoublePrecision xj, relgap(), delta1, delta2;
  
  Integer clustr_check();
  extern Integer count_list();
  extern Integer reduce_list2();
  extern Integer menode_ ();  
  extern Integer idamax_(), mxmynd_ (), mypid();
  extern void fil_int_lst();
  extern void fil_dbl_lst();
  
  /*
    intrinsic DoublePrecision precision
    */
  
  
  /*
    blas calls
    */
  
  extern void dcopy_(), dscal_(), daxpy_();
  extern DoublePrecision ddot_(), dasum_(), dnrm2_();
  
  /*
    lapack calls
    */
  
  extern DoublePrecision dlamch_();
  extern void xerbla_(), dlagts_(), dlagtf_(), dlarnv_();
  
  me = mxmynd_ ();
  num_eig = m;
  
  *num_clustr = 0;
  num_cls = I_ZERO;
  num_all_cls = I_ZERO;
  max_clustr_size = I_ZERO;
  c_ptr = clustr_info;

  
  /*
    assume that all eigenvalues are non-zero
    */
  
  /*
    this processor finds the eigenvectors between imin <= and <= imax
    */

  
  
  /*
    cluster information
    */
  
  for ( j = 0; j <  4 * n; j++ )
    c_ptr[j] = -1;
  
  c_ptr = clustr_info;
  
  /*
    should have the same error messages as lapack
    */
  
  /*
    get machine precision 
    */
  
  eps = DLAMCHE;
  

  *num_clustr = 1;
  b1 = 0;
  bn = n-1;
  *(c_ptr++) = b1;
  *(c_ptr++) = bn;
  *(c_ptr++) = b1;
  *(c_ptr++) = bn;
  return(bn-b1+1);
  
  /*
    
    should spread this across blocks;
    assume, of course, that nblk is set correctly
    
    */
  
  beg_of_block = -1;
  end_of_block = -1;
  
  while ( end_of_block < (num_eig-1 )) {
    ++end_of_block;
    beg_of_block = end_of_block;
    
    /*
       find the end of the block
       */
    
    
    
    while (( iblock[end_of_block] == iblock[beg_of_block] )
	   && ( end_of_block < num_eig )) {
      end_of_block++;
    }
    
    --end_of_block;
    
    /*
     */
    
    nblk = iblock[beg_of_block] - 1;
    
    if (nblk == 0 ) 
      b1 = 0 ;
    else
      b1 = isplit[nblk-1];
    
    bn = isplit[nblk]-1;
    blksiz = bn - b1 + 1;
    c1 = b1;
    delta1 = (DoublePrecision) max(100, n) ;
    delta2 = min(1.e-3, 1.e0/(DoublePrecision) n);
    
    /*
      Skip all the work if the block size is one.
      */
    
    if (blksiz == 1) {
      *(c_ptr++) = c1;
      *(c_ptr++) = c1;
      *(c_ptr++) = b1;
      *(c_ptr++) = bn;
      num_cls++;
    }

    c1 = b1;
    c2 = bn;
    
    if ( blksiz == 2 ) {
      if ( relgap(w[b1],w[bn]) < delta1 ) {
	*(c_ptr++) = c1;
	*(c_ptr++) = c2;
	*(c_ptr++) = b1;
	*(c_ptr++) = bn;
	num_cls++;
      }
      else {
	*(c_ptr++) = c1;
	*(c_ptr++) = c1;
	*(c_ptr++) = b1;
	*(c_ptr++) = bn;
	num_cls++;
	*(c_ptr++) = c2;
	*(c_ptr++) = c2;
	*(c_ptr++) = b1;
	*(c_ptr++) = bn;
	num_cls++;
      }
    }

    
    if ( blksiz > 2 ) {
      jblk = 0;
      ptr = beg_of_block;
      tail = beg_of_block;
      clustr_size = 1;
      xj = w[ptr];

      for (j = beg_of_block+1; j < end_of_block+1 ; ++j) {
	if ( j == (end_of_block )){
	  *(c_ptr++) = ptr;
	  *(c_ptr++) = j;
	  *(c_ptr++) = b1;
	  *(c_ptr) = bn;
	  num_cls++;
	  max_clustr_size = max(( ptr - tail + 1), max_clustr_size);
	  break;
	}
	
	if ( clustr_size == 1 ) {
	  xj1 = w[j];
	  if ( relgap(xj, xj1) < delta1 ) {
	    xj = w[j-1];
	    xj1 = w[j];
	    tail++;
	    clustr_size = 2;
	  }
	  else {
	    *(c_ptr++) = ptr;
	    *(c_ptr++) = tail;
	    *(c_ptr++) = b1;
	    *(c_ptr++) = bn;
	    num_cls++;
	    max_clustr_size = max(( ptr - tail + 1), max_clustr_size);
	    ptr = j;
	    tail = j;
	    xj = w[j];
	    clustr_size = 1;
	  }
	}
	else {
	  xj2 = w[j];
	  if ( min(relgap( xj, xj1),relgap(xj1, xj2)) > delta2 ) {
	    *(c_ptr++) = ptr;
	    *(c_ptr++) = tail;
	    *(c_ptr++) = b1;
	    *(c_ptr++) = bn;
	    num_cls++;
	    max_clustr_size = max( (ptr-tail+1) , max_clustr_size);
	    ptr = j;
	    tail = j;
	    xj = w[j];
	    clustr_size = 1;
	  }
	  else{
	    tail++;
	    clustr_size++;
	    xj = w[j-1];
	    xj1 = w[j];
	  }
	}
      }
    }
  }
  
  jjj =0;
  for ( iii = 0; iii < num_cls; iii++ ) {
    printf(" cptr c1 %d cn  %d b1 %d bn %d \n", clustr_info[jjj], clustr_info[jjj+1], clustr_info[jjj+2], clustr_info[jjj+3]);
    jjj += 4;
  }
  
  *num_clustr = num_cls;
  
  return(max_clustr_size);
  

  
}


