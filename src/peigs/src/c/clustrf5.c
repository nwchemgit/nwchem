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
 *********************************************
 
 PeIGS internal routine for eigen-value
 cluster analysis
 
      This version ASSUMES eigenvalues
      are accurate relative to norm(T), where T is the full
      tridiagonal matrix, but not necessarily relative to norm(Ti)
      if norm(Ti) < norm(T), where Ti is one of the submatrices into 
      which T breaks (as per isplit, nsplit).  If eigenvalues in block
      Ti are always accurate relative to norm(Ti), then ortol should be
      based on norm(Ti) rather than norm(T). 

      This assumption means that the orthogonalization criterion for
      all blocks is based on norm(T), not norm(Ti) for each block.  This
      means we do orthogonalization more often than we would if the 
      orthogonalization criterion for each block was based on norm(Ti).

      This assumption about accuracy is consistent with taking 
      abstol <= 0.0 in dstebz, pdspev, etc., i.e., then dstebz/pstebz
      computes eigenvalues with an error less than eps * norm(T).

      In principle the user could specify 0.0 < abstol < eps * norm(T) 
      in dstebz/pstebz
      in which case the assumption of this routines may be violated.
      Violation of this assumption will not reduce the accuracy of the
      results, it will only cause the code to run slower since it
      will cause more eigenvectors to be orthogonalized than is strictly
      necessary.

      If we do not make this assumption, and compute the orthogonalization
      criterion based on norm(Ti), but the eigenvalues are only accurate
      relative to norm(T) but not norm(Ti), then the computed eigenvectors
      are likely to have poor (possible no) orthogonality.
      
 
      Integer clustrf_ (n, d, e, m, w, mapZ, vecZ, iblock, nsplit, isplit,ptbeval, num_clustr, clustr_info, *imin,
 , proclist, iscratch )
 
 Integer *n, *m, *mapZ, *iblock, *nsplit, *isplit, *num_clustr, *clustr_info, *imin, *proclist, *iscratch;
 DoublePrecision *d, *e, *w, **vecZ, *ptbeval;
 
 
 not intended for separate call
 not protected from user input errors
 
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "globalp.c.h"

#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define MIN(a,b) ((a) < (b) ? (a) : (b))

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


Integer clustrf5_ (n, d, e, m, w, mapZ, vecZ, iblock, nsplit, isplit, ptbeval, num_clustr,
		   clustr_info,  imin, proclist , nacluster, icsplit, iscratch)
     Integer *n, *m, *mapZ, *iblock, *nsplit, *isplit, *num_clustr,
  *clustr_info, *imin, *proclist, *nacluster, *icsplit, *iscratch;
     DoublePrecision *d, *e, *w, **vecZ, *ptbeval;
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

       nacluster = total number of clusters in the list w of m
                     eigenvalues,  including clusters of size 1 in
                     blocks of size 1.

       icsplit ....... Length m.
                     icsplit[i] is 
                     the index (0 to m-1) of the last eigenvalue in the i-th
                     cluster in w, i = 0 to nacluster-1.
                     The i-th cluster, i = 0 to nacluster - 1 is
                     the set of eigenvalues

                     0 to icsplit[0],              i = 0
                     icsplit[i-1]+1 to icsplit[i], i = 1 to nacluster-1.
       
       */
{
  Integer jblk, nblk, msize;
  Integer i, j, num_cls, num_all_cls;
  Integer b1, num_eig;
  Integer max_clustr_size;
  Integer beg_of_block, end_of_block;
  Integer bn;
  Integer clustrptr=0, blksiz;
  Integer me;
  Integer *c_ptr;
  Integer iflag, k, imax;
  Integer ii, nn_proc;
  
  DoublePrecision tmp, *eval, sep, eps;
  DoublePrecision onenrm, pertol;
  DoublePrecision xjm, eps1;
  DoublePrecision ortol, xj, sepfine;
  
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
  
  extern void xerbla_(), dlagts_(), dlagtf_(), dlarnv_();
  
  me = mxmynd_ ();
  num_eig = *m;
  
  *num_clustr = 0;
  num_cls = 0;
  msize = *n;
  num_all_cls = I_ZERO;
  *nacluster = I_ZERO;
  max_clustr_size = I_ZERO;
  
  /*
      Compute reorthogonalization criterion.  

      ASSUMES eigenvalues
      are accurate relative to norm(T), where T is the full
      tridiagonal matrix, but not necessarily relative to norm(Ti)
      if norm(Ti) < norm(T), where Ti is one of the submatrices into 
      which T breaks (as per isplit, nsplit).  If eigenvalues in block
      Ti are always accurate relative to norm(Ti), then ortol should be
      based on norm(Ti) rather than norm(T) (the commented out computation
      of ortol latter in the code was used to do this before even though
      the eigenvalues weren't really accurate relative to norm(Ti). ).
   */
      
  onenrm = fabs( d[0] ) + fabs( e[1] );
  tmp = fabs(d[msize-1]) + fabs(e[msize-1]);
  onenrm = MAX(onenrm, tmp);
  for (i = 1; i < msize-1; ++i) {
    tmp = fabs(d[i]) + fabs(e[i]) + fabs(e[i + 1]);
    onenrm = MAX(onenrm, tmp);
  }
  
  ortol = onenrm * (DoublePrecision ) 1.e-3 ;
  sepfine = MAX(1.e2, (DoublePrecision) msize) * sqrt(DLAMCHE);
  sepfine = MIN(sepfine, MIN(1.e-3, 1.0e0/((DoublePrecision) msize)));
  
  /*
    sepfine = sepfine*MAX(ortol, 1.e0);
    */
  
  c_ptr = iscratch;
  nn_proc = reduce_list2( num_eig, mapZ, c_ptr);
  
  iflag = -1;
  for ( j = 0; j < nn_proc; j++ ){
    if ( c_ptr[j] == me ) {
      iflag = j;
      break;
    }
  }    
  
  k = 0;
  imax = -1;
  *imin = -1;
  for ( j = 0; j < iflag; j++ ) {
    imax = count_list(c_ptr[j], mapZ, &num_eig);
    k += imax;
  }	

  
  /*
    this processor finds the eigenvectors between imin <= and <= imax
    */

  if ( iflag == -1 ){
    *imin = -1;
    imax = -1;
  }
  else {
    *imin = k;
    imax = k + count_list( me, mapZ, &num_eig) - 1;
  }
  
  c_ptr = iscratch;
  nn_proc = reduce_list2(num_eig, mapZ, c_ptr );
  
  c_ptr += nn_proc;
  for ( j = 0; j < num_eig; j++ )
    *(c_ptr++ ) = mapZ[j];
  
  c_ptr = iscratch + nn_proc;
  
  k = 0;
  for ( j = 0; j < nn_proc; j++ ) {
    ii = iscratch[j];
    end_of_block = count_list(ii, c_ptr, &num_eig);
    fil_int_lst (end_of_block, &mapZ[k], ii);
    k += end_of_block;
  }
  
  eval = ptbeval;
  c_ptr = clustr_info;
  
  /*
    cluster information
    */
  
  for ( j = 0; j <  4 * msize; j++ )
    *(c_ptr++) = -1;
  
  for ( j = 0; j < num_eig; j++ )
    icsplit[j] = -5;

  if ( iflag == -1 )
    return(-50);
  
  c_ptr = clustr_info;
  
  /*
    should have the same error messages as lapack
    */
  
  /*
    get machine precision 
    */
  
  eps = DLAMCHE;    
  
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
    
    
    
    while (( iblock[end_of_block] == iblock[beg_of_block] ) && ( end_of_block < num_eig )) {
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
    
    /*
      Skip all the work if the block size is one.
      */
    
#ifdef DEBUG1
    fprintf(stderr, "me = %d b1 = %d  bn = %d blksiz1 = %d \n", me, b1, bn, blksiz );    
#endif
    
    if (blksiz == 1) {
      if ( clustr_check(beg_of_block, end_of_block, *imin, imax) == 1 ) { 
	ii = beg_of_block - *imin;
	fil_dbl_lst ( msize, vecZ[ii], 0.0e0);
	vecZ[ii][b1] = 1.0e0;
      }
      
      icsplit[ num_all_cls ] = end_of_block;
      num_all_cls++;
      *nacluster = num_all_cls;
      
    }
    
#ifdef DEBUG1
    fprintf(stderr, " got here 1 me = %d \n", me );
#endif

    if ( blksiz > 1 ) {
      /*
	Compute reorthogonalization criterion and stopping criterion
	Gershgorin radius
	
	should use pointers or loop unrooling here to speed things up; this is a O(n) operation
	so things aren't that expensive
	
	*/
      
/*
 *  This is how one would compute ortol if the eigenvalues of block i
 *  of T, Ti, were accurate relative to norm(Ti).
 *
 *     onenrm = fabs( d[b1] ) + fabs( e[b1+1] );
 *     tmp = fabs(d[bn]) + fabs(e[bn]);
 *     onenrm = max(onenrm, tmp);
 *     for (i = b1 + 1; i < bn; ++i) {
*	tmp = fabs(d[i]) + fabs(e[i]) + fabs(e[i + 1]);
*	onenrm = max(onenrm, tmp);
*      }
 *     
 *     ortol = onenrm * (DoublePrecision ) 1.e-3 ;
 */
      
#ifdef DEBUG1
    fprintf(stderr, " got here 2 me = %d \n", me );
#endif

      
      /*
	Loop through eigenvalues of block nblk.
	*/
      
      jblk = 0;
      for (j = beg_of_block; j < ( end_of_block + 1) ; ++j) {
	++jblk;
	xj = w[j];
	
	/*
	  If eigenvalues j and j-1 are too close, add a relatively small
	  perturbation.  If the eigenvalues are not in increasing 
	  order in the block, exit.
	*/
	
	if (jblk > 1) {
	  eps1 = fabs(eps * xj);
	  pertol = eps1 * R_TEN;
	  if ((xj - w[j-1]) < -eps1) {
	    printf(" Error in ordering eigenvalues: -5 error clustrf me = %d \n",
                   (int) me );
	  }
	  sep = xj - xjm;
/*
	  if (sep < pertol*MAX(fabs(xj), fabs(xjm)))
*/
	if ( sep < pertol)
	 xj = xjm + pertol;
	}
	eval[j] = xj;
	
#ifdef DEBUG1
	fprintf(stderr, " got here 3 me = %d \n", me );
#endif
	
	if (jblk == 1)
	  clustrptr = j;
	
	/*
	  tight cluster
	  */
	
	if ( jblk > 1 )  {            /* jblk > 1 */
	  sep = fabs(xj - xjm);
	  /*
	    printf(" got here 4 me = %d sep xj %20.16g xjm %20.16g  sep %20.16g \n", me, xj, xjm, sep-sepfine*MAX(fabs(xj), fabs(xjm)));
	  */
	  
	  if (fabs(sep) > fabs(sepfine*MAX(fabs(xj),fabs(xjm)))) {
	    if ( clustr_check(clustrptr, j-1, *imin, imax) == 1 ) {
	      *(c_ptr++) = clustrptr;
	      *(c_ptr++) = j-1;
	      *(c_ptr++) = b1;
	      *(c_ptr++) = bn;
	      num_cls++;
	      max_clustr_size = MAX( j - clustrptr, max_clustr_size);
	    }
	    clustrptr = j;
	    jblk = 1;

            icsplit[ num_all_cls ] = j-1;
            num_all_cls++;
            *nacluster = num_all_cls;


#ifdef DEBUG1
	    fprintf(stderr, " out of 4 me = %d \n", me );
#endif

	  }
	}
	
	/* 
	  assume that xj - xjm < ortol but we're at the end of the blk 
	  */
	if ( j == end_of_block ) {
	  if ( clustr_check(clustrptr, end_of_block, *imin, imax) == 1 ) {
	    
#ifdef DEBUG1
	    fprintf(stderr, " got here 5 me = %d \n", me );
#endif
	    
	    *(c_ptr++) = clustrptr;
	    *(c_ptr++) = end_of_block;
	    *(c_ptr++) = b1;
	    *(c_ptr++) = bn;
	    num_cls++;
	    max_clustr_size = MAX( (end_of_block -clustrptr + 1), max_clustr_size);
#ifdef DEBUG1
	    fprintf(stderr, " out of 5 me = %d cbeg %d cend %d b1 %d bn %d \n",me,
		    clustrptr, end_of_block, b1, bn );
#endif
	  }

          icsplit[ num_all_cls ] = end_of_block;
          num_all_cls++;
          *nacluster = num_all_cls;

	}
	
	xjm = xj;
      }
    }
  }
  
  /*
  ime = -1;
  ii = -1;
  for ( ii = 0; ii < msize; ii++ ){
    if ( mapZ[ii] == me ){
      ime = ii;
      break;
    }
  }

  if ( ime % 2 == 1 ) {
    if ( num_all_cls > 1 ) {
      iscratch[0] = clustr_info[0];
      iscratch[1] = clustr_info[1];
      iscratch[2] = clustr_info[2];
      iscratch[3] = clustr_info[3];
      c_ptr = &clustr_info[ 4 * ( num_all_cls - 1)];
      clustr_info[0] = c_ptr[0];
      clustr_info[1] = c_ptr[1];
      clustr_info[2] = c_ptr[2];
      clustr_info[3] = c_ptr[3];
      c_ptr = &clustr_info[0];
      c_ptr[0] = iscratch[0];
      c_ptr[1] = iscratch[1];
      c_ptr[2] = iscratch[2];
      c_ptr[3] = iscratch[3];
    }
  }
  */
  
  c_ptr = clustr_info;
  *num_clustr = num_cls;

  
/*
  for ( ii = 0; ii < 4* *num_clustr; ii++ )
    printf("me = %d clustr_info[%d] =  %d \n", me, ii, *(c_ptr++));
*/
  
  /*
    for ( ii = 0; ii < num_eig ; ii++ )
    printf("me = %d mapZ[%d] =  %d \n", me, ii, mapZ[ii]);
    
    printf("me = %d exiting clustrf2 imin =  %d \n", me, *imin);
    */
  
  
/*
    for ( ii = 0; ii < 4* *num_clustr; ii++ )
    printf("me = %d clustrf5 close clustr_info[%d] =  %d \n", me, ii, *(c_ptr++));
    
*/
  
    return(max_clustr_size);
}

