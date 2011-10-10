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
/*
  PeIGS internal routine:
  
*/

#include <stdio.h>
#include <math.h>

#include "globalp.c.h"
#include "clustr_inv.h"

/*
  #define DEBUG1
  #define DEBUG5
  */
#define CLUSTRLEN  4
#define LOOP  3
#define INV_TIME 2
#define ITIME1  2

Integer clustrinv4_(n, d, e, dplus, lplus, ld, lld, eval, schedule, num_clustr, mapZ, mapvecZ, vecZ, imin, nacluster, icsplit, iscratch, scratch)
     Integer *n, *schedule, *num_clustr, *mapZ, *mapvecZ, *imin,
  *nacluster, *icsplit, *iscratch;
     DoublePrecision *d, *e, *ld, *lld, *eval, **vecZ, *scratch, *dplus, *lplus;
{
  /*
    n = dimension of the tridiagonal matrix
    d = diagonal of the tridiagonal matrix
    e[2:n] = upper diagonal of the tridiagonal matrix
    schedule = a integer array of the scheduling information
    mapZ = location of the eigenvectors
    mapvecZ = in
    vecZ
    iscratch = integer scratch space
    scratch = DoublePrecision precision scratch space
    
    returns the number of eigenvector that this processor holds
    */
  
  /*
  static Integer three = 3, IONE = 1;
  */

  Integer bb1, bn, iii, ii;
  Integer clustr_ptr, cn;
  Integer me, naproc, Zvec;
  Integer *cl_ptr;
  Integer c1, csiz;
  Integer cl_num;
  Integer ibad;
  
  DoublePrecision eps;
  DoublePrecision *dscrat;
  
  extern void xerbla_();
  extern Integer idamax_(), mclock_(), succ_(), mxmynd_(), mxnprc_();
  extern DoublePrecision dasum_(), dnrm2_(), ddot_(), dlarnd_();
  extern void dlagtf_(), dlarnv_();
  extern void printff_(), mgs_3();
  extern void mgspnl_(), dscal_(), dlagts_();
  extern void dcopy_(), daxpy_();
  extern Integer count_list ();
  extern void mgs (), mgscs();
  extern void fil_dbl_lst ();
  
  /*
    Function Body
    */
  
  me = mxmynd_();
  naproc = mxnprc_();

  ibad = 0;

  dscrat = scratch;

  /*
    Get machine constants. should set this up somewhere to call it only once
    */

  eps = DLAMCHE;
  
  /*
    Initialize seed for random number generator DLARND.
    how about randomizing the randomness?
    */
  
  cl_num = *num_clustr;
  Zvec = 0;
  
  if ( cl_num == 0 )
    return(ibad);
  
  /*
   * If I own the start of my last cluster, but not all of this cluster
   * then do my part of the last cluster and send it to the processor
   * which owns the first vector in the cluster which is not mine.
   */
  
  cl_ptr = schedule;
  ibad = 0;

  /*
    bb1 = 0;
    bn = *n-1;
    c1 = 0;
    cn = *n-1;
    mgscs( n, vecZ, mapZ, bb1, bn, c1, cn, iscratch, dscrat);
    return(ibad);
    */
  
  
  for (clustr_ptr= 0;  clustr_ptr < cl_num ; clustr_ptr++) {
    c1 = *(cl_ptr++);
    cn = *(cl_ptr++);
    bb1 = *(cl_ptr++);
    bn = *(cl_ptr++);
    csiz = cn -c1 +1;
    
    iii = 0;
    for ( ii = c1; ii <= cn; ii++ ){
      if ( mapZ[ii] == me ) {
	iii = 1;
	break;
      }   
    }
    
    if ((iii == 1 ) && (csiz > 1 )) {
      mgscs( n, vecZ, mapZ, bb1, bn, c1, cn, iscratch, dscrat);
    }
  }
  
  /*
    printf( " me = %d Exiting clustr4_ \n", me );
    fflush(stdout);
    */
  
  return(ibad);
}

/* $Id$ */
