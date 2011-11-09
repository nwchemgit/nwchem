/*
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
/* **********************************************
   
   PeIGS internal routine
   
   Integer inv_it4( iii, n, c1, cn, b1, bn, Zbegin, map, mapvec, vector, d, e, eval, eps, stpcrt, onenrm, iwork, work)
   Integer iii, *n, *c1, *cn, *b1, *bn, *Zbegin, *map, *mapvec, *iwork;
   DoublePrecision *d, *e, **vector, *eval, *work, *eps, *stpcrt, *onenrm;
   
   This routine performs inverse iteration.
   
   */   

#include <stdio.h>
#include <memory.h>
#include <math.h>

#include "globalp.c.h"

#include "clustr_inv.h"

#define    ITS   1
#define    MAXITR 100

Integer inv_it4(iii,  n, c1, cn, b1, bn, Zbegin, map, mapvec, vector, d, e, ld, lld, eval, eps, stpcrt, onenrm, iwork, work)
     Integer *iii, *n, *c1, *cn, *b1, *bn, *Zbegin, *map, *mapvec, *iwork;
     DoublePrecision *d, *e, *ld, *lld, **vector, *eval, *work, *eps, *stpcrt, *onenrm;
     
     /*
       this routine performs inverse iteration on *n vectors given by map[0:*n-1]
       and whose storage location on each processor is given by mm
       
       n = dimension of the vectors
       m = number of vectors to inverse iterate
       
       d is the diag of the matrix
       e is the sub and super diagonals of the matrix
       eval is the list of the eigenvalues  eval(c1) ... eval[cn] are the eigenvalues; mapvec[k] gives
       the true index of the eigenvalues
       
       */
{
  static Integer IONE = 1;
  Integer j, blksz, i_1;
  Integer indrv1, indrv2, indrv3, indrv4, indrv5;
  Integer k, lsiz, me;
  DoublePrecision xj, sep;
  Integer mxmynd();
  Integer csiz, ibad;
  Integer i;
  DoublePrecision delta = 0.;
  DoublePrecision ztz, *wwork1;
  Integer zbegin1, kk, zend, bb1 = *b1+1, bbn = *bn+1, *iwork1;
  
  extern void dscal_(), dlagtf_(), dlagts_(), dcopy_(), dscal_();
  extern Integer idamax_();
  extern DoublePrecision dnrm2_(), dasum_(), ddot_();
  extern Integer mxmynd_();
  extern void dgetavec2_();
  extern void dgetavec3_();
  
  extern void mgs_prev ();
  
  me = mxmynd_();
  
  blksz = *bn - *b1 + 1;
  
  lsiz = *n;
  /*
  iwork1 = (Integer *) malloc( 16* lsiz *sizeof(Integer));
  wwork1 = (DoublePrecision *) malloc( 16* lsiz *sizeof(DoublePrecision));
  */
  iwork1 = iwork;
  wwork1 = work;
  
  for(i = 0;i < 16*lsiz;i++){
    iwork1[i] = 0;
    wwork1[i] = 0.;
  }
  
  indrv1 = 0;
  indrv2 = lsiz + indrv1;
  indrv3 = lsiz + indrv2;
  indrv4 = lsiz + indrv3;
  indrv5 = lsiz + indrv4;
  
  ibad = 0;
  csiz = *cn - *c1 + 1;
  
#ifdef DEBUG  
  fprintf(stderr, " csiz = %d blksize = %d \n", csiz, blksz);
#endif
  
  /*
    iwork1 = iwork+*n;
    */
  
  sep = 0.;
  k = *Zbegin;  
  i_1 = *c1;
  
  for ( j = 0; j < csiz ;  j++ ) {
    if ( map[ i_1 ] == me ) {
      mapvec[k] = i_1;
      xj = eval[i_1];
      if ( *iii == 0 ){
    	dgetavec2_( &j, &xj, &delta, n, &bb1, &bbn, e, d, ld, lld, wwork1,
		    wwork1 + lsiz, wwork1 + 2*lsiz, wwork1 + 3*lsiz, wwork1 + 4*lsiz, wwork1 +
		    5* lsiz, wwork1 + 6* lsiz, vector[k], &kk, &ztz, &zbegin1, &zend, &j, iwork1 );
	
      }
      else {
	delta = fabs(xj - eval[i_1-1] )/fabs(xj) * 2.e-16;
	dgetavec3_( &xj, &delta, n, &bb1, &bbn, e, d, ld, lld, wwork1,
		    wwork1 + lsiz, wwork1 + 2*lsiz, wwork1 + 3*lsiz, wwork1 + 4*lsiz, wwork1 +
		    5* lsiz, wwork1 + 6* lsiz, vector[k], &kk, &ztz, &zbegin1, &zend, &j, iwork1 );
      }
      ztz = 1./ztz;      
      /*
	printf(" inv_it4 me = %d xj = %f j = %d c1 = %d cn = %d ztz %f k = %d \n", me,  xj, j, *c1, *cn, ztz, k);
      */
      dscal_(&blksz, &ztz, vector[k] + *b1, &IONE);
      i_1++;
      k++;
    }
  }
  
  /*
  free(wwork1);
  free(iwork1);
  */
  printf(" exit inv_it4 me = %d \n", (int)me);  
  return(ibad);
}  


/* $Id$ */
