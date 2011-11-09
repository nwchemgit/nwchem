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

/* **********************************************

   PeIGS internal routine

Integer inv_it( n, c1, cn, b1, bn, Zbegin, map, mapvec, vector, d, e, eval, eps, stpcrt, onenrm, iwork, work)
     Integer *n, *c1, *cn, *b1, *bn, *Zbegin, *map, *mapvec, *iwork;
     DoublePrecision *d, *e, **vector, *eval, *work, *eps, *stpcrt, *onenrm;

     This routine performs inverse iteration.
     
     */   




#include <stdio.h>
#include <math.h>

#include "globalp.c.h"

#include "clustr_inv.h"

#define    ITS   2
#define    MAXITR 10
#define FABS(a)   ((a) < (DoublePrecision) 0.0e0 ? (- a): (a))

extern DoublePrecision psigma, psgn;


Integer inv_it3( n, c1, cn, b1, bn, Zbegin, map, mapvec, vector, d, e, eval, eps, stpcrt, onenrm, iwork, work)
     Integer *n, *c1, *cn, *b1, *bn, *Zbegin, map[], mapvec[], iwork[];
     DoublePrecision d[], e[], **vector, eval[], work[], *eps, *stpcrt, *onenrm;
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
  static Integer IONE = 1, IMINUSONE = -1;
  Integer j, blksz, i_1, niter;
  Integer nrmchk, indrv1, indrv2, indrv3, indrv4, indrv5;
  Integer jmax, k, l, me;
  DoublePrecision xj, nrm, scl;
  DoublePrecision tol, *ptr;
  Integer info, mxmynd();
  Integer csiz, ibad;
  Integer indx22;
  
  extern void dscal_(), dlagtf_(), dlagts_(), dcopy_(), dscal_();
  extern Integer idamax_();
  extern DoublePrecision dnrm2_(), dasum_(), ddot_();
  extern Integer mxmynd_();
  
  extern void mgs_prev ();
  
  me = mxmynd_();
  
  blksz = *bn - *b1 + 1;
  
  l = *n;
  
  indrv1 = 0;
  indrv2 = l + indrv1;
  indrv3 = l + indrv2;
  indrv4 = l + indrv3;
  indrv5 = l + indrv4;
  
  ibad = 0;

  l = 0;
  csiz = *cn - *c1 + 1;
  k = *Zbegin;
  
#ifdef DEBUG  
  fprintf(stderr, " csiz = %d k = %d \n", csiz, k);
#endif
  
  
  /*
    Compute LU factors with partial pivoting  ( PT = LU )
  */
  
  
  for ( j = 0; j < csiz ;  j++ ) {
    i_1 = j+*c1;
    if ( map[ i_1 ] == me ) {
      mapvec[k] = i_1;
      xj = eval[i_1];
      xj += psgn*psigma;
      
#ifdef DEBUG  
      fprintf(stderr, "xj = %g \n", xj);
#endif      
      
      ptr = &work[indrv1];
      dcopy_(&blksz, &vector[k][*b1], &IONE, ptr, &IONE );
      dcopy_(&blksz, &d[*b1], &IONE, &work[indrv4], &IONE);
      i_1 = blksz - 1;
      dcopy_(&i_1, &e[*b1+1], &IONE, &work[ indrv2 + 1 ], &IONE);
      dcopy_(&i_1, &e[*b1+1], &IONE, &work[ indrv3 ], &IONE);
      
      tol = ZERO;
      dlagtf_( &blksz, &work[indrv4], &xj, &work[indrv2 + 1], &work[indrv3], &tol,
	       &work[indrv5], &iwork[0], &info);
      
      /*
        Update iteration count.
	*/
      
      nrmchk = 0;
      nrm = ZERO;
      
      /*
	loop for inverse iteration
	*/
      
      niter = 0;
      
      loop :;
      
      niter++;
      
      if (niter > MAXITR){
	fprintf(stderr, " \n me = %d MAXITS exceeded in inv_it \n \n", (int)me);
	fprintf(stderr, " c1= %d cn= %d b1= %d bn = %d \n",
                (int)*c1,(int)*cn,(int)*b1,(int)*bn);
	fprintf(stderr, " nrm = %g stpcrt = %g \n", nrm, *stpcrt );
	fprintf(stderr, "xj = %e onenrm = %e \n", xj, *onenrm);
	for (indx22  = *b1; indx22 < *bn + 1; indx22++ )
	  fprintf(stderr, "d[%d] = %e e[%d] = %e \n",
		  (int)indx22, d[indx22], (int)indx22, e[indx22]);

	
	printf(" \n me = %d Eigenvector %d failed to convergve in inv_it3 \n \n",
               (int)me, (int)(*c1 + j) );
	
        if( ibad == 0 )
	  ibad = *c1 + j + 1;   /* note +1 needed since *c1 may be 0 */
	
        continue;
      }
      
      
      /*
	Normalize and scale the righthand side vector Pb.
      */
      
      scl = dasum_(&blksz, ptr, &IONE);
      scl = ((DoublePrecision) 1.0e0)/scl;
      scl *= (DoublePrecision) blksz ;
      scl *= *onenrm ;
      /*
	scl *= max(*eps, FABS(work[indrv4 + blksz - 1]));
      */

      scl *= max(*eps, FABS(work[indrv4 + blksz - 1]) / *onenrm );
      dscal_(&blksz, &scl, ptr, &IONE );
      
      /*
	Solve the system LU = Pb.
	*/
      
      tol = ZERO;
      
      
      dlagts_( &IMINUSONE, &blksz, &work[indrv4], &work[ indrv2 + 1],
	       &work[indrv3], &work[indrv5], &iwork[0], ptr, &tol, &info);
      
      jmax = idamax_(&blksz, ptr, &IONE);
      nrm = FABS(ptr[jmax - 1]);
      
      if ( nrm < *stpcrt ) goto loop;
      nrmchk++;

      if ( nrmchk < ITS ) goto loop;
      /*
	should do something about this

      if ( nrmchk > 5 ) 
	exit(1);
	*/
      
      scl = dnrm2_(&blksz, ptr, &IONE);
      scl = ((DoublePrecision) 1.0e0)/scl;
      dscal_(&blksz, &scl, ptr, &IONE);
      dcopy_( &blksz, ptr, &IONE, vector[k] + *b1, &IONE);
      k++;
    }
  }
  return(ibad);
}  


