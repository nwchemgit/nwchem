/*======================================================================
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
  
  not intended for external use; see the routine pstein
  
  Integer clustrinv_(n, d, e, eval, schedule, num_clustr, mapZ, mapvecZ, vecZ, imin, nacluster, icsplit, iscratch, scratch)
  
  Integer *n, *schedule, *num_clustr, *mapZ, *mapvecZ, *Zbegn, *nacluster,
  *icsplit, *iscratch;
  
  DoublePrecision *d, *e, *eval, **vecZ, *scratch;
  
  
  for computing eigenvectors by inverse iteration with cluster modified Gram Schmidt
  with "de-tangling" of clusters.

  returns 0     if everything seems ok.
          k > 0 if k-th eigenvector (owned by processor mapZ[k-1], k = 1 to neigval,
                failed to converge in inverse iteration.  In this case eigenvectors
                1 to k-1 converged ok, eigenvector k failed to converge, and it is
                possible that some eigenvectors numbered greater than k failed to
                converge.  In any case, the first k-1 eigenvectors are good,
                while eigenvectors k is bad, and eigenvectors k+1 to neigval are
                probably bad too.  Note that the meaning of bad is not always clear,
                in general, however, it means that the eigenvector is probably
                not accurate to full precision.
                

  This code "de-tangles" clusters as follows.  Assume the eigenvalues in a cluster
  are owned by more than one processor (as determined by mapZ).  Then, the processor
  owing the first eigenvector in the cluster does all required iterations of
  inverse iteration and mgs to compute the eigenvectors in this cluster.  It then
  sends these finished eigenvectors to the second processor in mapZ which owns part
  of this cluster.  At this point the clusters are de-tangled and all processors can
  start working on their eigenvectors.
  
  This version of the code assumes that during de-tangling each processor owns the
  first eigenvector in at most one multi-processor cluster  and receives the finished
  eigenvalues from at most one other processor.  This assumptions are satified by
  the current algorithm in clustrf and are checkd by this routine.

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
  
  static Integer three = 3, IONE = 1;
  Integer indx, i, j, iseed[4], bb1, bn, jjj, iii, ii;
  Integer blksiz, clustr_ptr, cn;
  Integer me, naproc, Zvec;
  Integer *cl_ptr;
  Integer c1, csiz, xc1, xcsiz, xblksiz;
  Integer cl_num;
  Integer itime;
  Integer send_num, send_cl, send_to,
          recv_num, recv_cl, recv_from,
          myindx, ime, itype, nvecs, isize, ival, first, ibad, itmp;
  
  DoublePrecision stpcrt, onenrm, eps;
  DoublePrecision tmp, *dscrat, *first_buf;
  
  extern void xerbla_();
  extern Integer idamax_(), mclock_(), succ_(), mxmynd_(), mxnprc_();
  extern DoublePrecision dasum_(), dnrm2_(), ddot_(), dlarnd_();
  extern void dlagtf_(), dlarnv_();
  extern void printff_(), mgs_3();
  extern void mgspnl_(), dscal_(), dlagts_();
  extern void dcopy_(), daxpy_();
  extern Integer count_list ();
  extern Integer inv_it();
  extern void mgs ();
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

  cl_ptr = schedule;
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
