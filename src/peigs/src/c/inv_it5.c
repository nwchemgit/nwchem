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
   
   This routine performs inverse iteration for bidiagonal matrices.
   
   */   

#include <stdio.h>
#include <math.h>

#include "globalp.c.h"

#include "clustr_inv.h"

#define    ITS   1
#define    MAXITR 100

Integer inv_it5( n, c1, cn, b1, bn, Zbegin, map, mapvec, vector, d, e, ld, lld, eval, eps, stpcrt, onenrm, iwork, work)
     Integer *n, *c1, *cn, *b1, *bn, *Zbegin, map[], mapvec[], iwork[];
     DoublePrecision d[], e[], ld[], lld[], **vector, eval[], work[], *eps, *stpcrt, *onenrm;
     
     /*
       this routine performs inverse iteration on *n vectors given by map[0:*n-1]
       and whose storage location on each processor is given by mm
       
       n = dimension of the vectors
       m = number of vectors to inverse iterate
       
       d is the diag of the bidiagonal matrix
       e is the super the bidiagonal diagonals of the matrix
       eval is the list of the eigenvalues  eval(c1) ... eval[cn] are the eigenvalues; mapvec[k] gives
       the true index of the eigenvalues
       
       */
{
  static Integer IONE = 1;
  Integer j, blksz, i_1;
  Integer k, lsiz, me;
  DoublePrecision xj, sep;
  Integer mxmynd();
  Integer csiz, ibad;
  DoublePrecision delta = 0., *p, *gamma, *dminus, *lplus, *t,
    *uminus, *dplus;
  DoublePrecision ztz=0.;
  Integer zbegin1, kk, zend, bb1 = *b1+1, bbn = *bn+1, *iwork1;
  
  extern void dscal_(), dlagtf_(), dlagts_(), dcopy_(), dscal_();
  extern Integer idamax_();
  extern DoublePrecision dnrm2_(), dasum_(), ddot_();
  extern Integer mxmynd_();
  extern DoublePrecision dgetavec_();
  
  
  extern void mgs_prev ();
  
  me = mxmynd_();
  
  blksz = *bn - *b1 + 1;
  
  lsiz = blksz;
  csiz = *cn - *c1 + 1;
  
  /*
    iwork1 = (Integer *) malloc( 16* lsiz *sizeof(Integer));
    lplus = (DoublePrecision *) malloc( 16*lsiz *sizeof(DoublePrecision));
  */
  
  iwork1 = iwork;
  lplus = work;

  dplus = lplus + lsiz;
  uminus = dplus + lsiz;
  dminus = uminus + lsiz;
  t = dminus + lsiz;
  p = t + lsiz;
  gamma = p + lsiz;
  
  ibad = 0;
  

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
      dgetavec_( &j, &xj, &delta, n, &bb1, &bbn, e, d, ld, lld, lplus,
		 dplus, uminus, dminus, t, p, gamma,
		 vector[k], &kk, &ztz, &zbegin1, &zend, &j, iwork1 );
      
#ifdef DEBUG
      printf(" invit5 me = %d csiz = %d j = %d xj = %f c1 = %d cn = %d ztz = %f \n", me, csiz, j, xj, *c1, *cn, ztz);
#endif
      
      ztz = 1.e0/ztz;
      dscal_(&blksz, &ztz, &vector[k][*b1], &IONE);
      i_1++;
      k++;
    }
  }
  
  /*
    free(lplus);
    free(iwork1);
  */
  
  return(ibad);
}  


/* $Id$ */
