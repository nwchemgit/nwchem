/*
 * $Id$
 *
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
#include <stdlib.h>

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

Integer clustrinv5_(n, d, e, dplus, lplus, ld, lld,
		    eval, r_eval, schedule, num_clustr, mapZ,
		    mapvecZ, vecZ, imin, nacluster,
		    icsplit, iscratch, scratch)
     Integer *n, *schedule, *num_clustr, mapZ[],
  *mapvecZ, *imin, *nacluster, *icsplit, *iscratch;
     DoublePrecision d[], e[], ld[], lld[], eval[], r_eval[], **vecZ,
  scratch[], dplus[], lplus[];
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
  
  static Integer IONE = 1;
  Integer indx, i, j, iseed[4], bb1, bn;
  Integer jjj;
  Integer blksiz, clustr_ptr, cn;
  Integer me, naproc, Zvec;
  Integer *cl_ptr;
  Integer c1, csiz, xc1, xcsiz, xblksiz;
  Integer cl_num;
  Integer itime, msize;
  Integer send_num, send_cl, send_to,
          recv_num, recv_cl, recv_from,
          myindx, ime, itype, nvecs, isize, ival, first, ibad, itmp;
  
  DoublePrecision stpcrt, onenrm, eps;
  DoublePrecision tmp, *dscrat, *first_buf;
  
  
  extern Integer mclock_(), succ_(), mxmynd_(), mxnprc_();
  extern void printff_();
  extern void mgspnl_();
  extern Integer inverm_();
  
  /*
    Function Body
    */
  
  me = mxmynd_();
  naproc = mxnprc_();
  msize = *n;

  first_buf = NULL;
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
  
  blksiz = 1;
  for (i = 0; i < 4; ++i)
    iseed[i] = 1;
  
  iseed[3] = 1;
  
  /*
    Compute eigenvectors of matrix blocks.
    */
  
  /*
    this is the number of clusters that this processor will process
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

  send_num = 0; 
  recv_num = 0;
  recv_cl = 0;
  send_cl = 0;
  send_to = -100;
  recv_from = -100;
  
  if( naproc > 1 ) {
    for (clustr_ptr= 0;  clustr_ptr < cl_num ; clustr_ptr++) {

      c1 = *(cl_ptr++);
      cn = *(cl_ptr++);
      bb1 = *(cl_ptr++);
      bn = *(cl_ptr++);
      
      if( cn > c1 ) {
        if( mapZ[c1] == me ) {
          for (j = c1+1; j <= cn; j++ ){
            if ( mapZ[j] != me ) {
              send_num++;
              send_cl = clustr_ptr;
              send_to = mapZ[j];
              break;
            }
          }
        }
        else {
        
          for (j = c1+1; j <= cn; j++ ){
            if ( mapZ[j] != mapZ[c1] ) {
              if ( mapZ[j] == me ) {
                recv_num++;
                recv_cl   = clustr_ptr;
                recv_from = mapZ[c1];
              }
              break;
            }
          }
        }
      }
    }

    /*
     *
     * The following tests verify that the "de-tangling" assumptions about
     * sending and/or receiving at most one finished set of eigenvectors listed
     * in the subroutine header are satisfied.
     *
     * The exit(-1) conditions should only occur if the
     * cluster scheduling algorithm is changed from its current form without
     * modifying the de-tangeling code appropriately.  Thus, the following exits
     * should never be executed in the tested, distribution version of the code.
     * Thus, we do not try to do a global exit here.
     *
     */

    if( send_num > 1 || recv_num > 1 ) {
      fprintf( stderr, " me = %ld Internal Error in PEIGS clustrinv. \n", me );
      fprintf( stderr, " me = %ld recv_num and/or send_num > 1. \n", me);
      exit( -1 );
    }
    if( ( send_num > 0 ) && ( send_cl != cl_num - 1) ) {
        fprintf( stderr, " me = %ld Internal Error in PEIGS clustrinv. \n", me );
        fprintf( stderr, " me = %ld send_cl != cl_num - 1. \n", me);
        exit( -1 );
    }
    if ( recv_num > 0 ){
	if ( recv_cl != 0 ) {
        fprintf( stderr, " me = %ld Internal Error in PEIGS clustrinv. \n", me );
        fprintf( stderr, " me = %ld recv_cl != 0. \n", me);
        exit( -1 );}
    }
  }

  if( recv_num > 0 || send_num > 0 ) {

    /* Set iscratch to a linked-list such that
     * iscratch[j] = k means processor
     * j receives from processor k when sending
     * the first block of a cluster to the next processor in that cluster.
     * k = -1 means processor j receives from no one.
     * Only really care about the part of the linked-list
     * including 'me'.
     */

    for (j= 0;  j < naproc; j++)
      iscratch[j] = -1;

    cn = -1;
    for( i = 0; i < *nacluster; i++ ) {
      c1 = cn + 1;
      cn = icsplit[i];
      if( cn > c1 ) {
        for (j = c1+1; j <= cn; j++ )
          if ( mapZ[j] != mapZ[c1] ) {
            iscratch[ mapZ[j] ] = mapZ[c1];
            break;
          }
      }
    }

    myindx = 0;
    j      = me;
    for (i= 0;  i < naproc; i++) {
      if( iscratch[j] == -1 ) 
        break;

      myindx++;
      j = iscratch[j];
    }

    ime = -1;
    for (i= c1;  i <= cn; i++) {
      ime ++;
      if ( mapZ[i] == me )
	break;
    }

    /*
      printf("me = %d send_to %d send_num %d recv_from %d recv_num %d \n", me, send_to, send_num, recv_from, recv_num);
      fflush(stdout);
      */
    

    if( iscratch[j] != -1 ) {
      fprintf( stderr, " me = %ld Internal Error in PEIGS clustrinv. \n", me );
      fprintf( stderr, " me = %ld Swapping of initial cluster data \n", me);
      fprintf( stderr, " me = %ld does not have a well defined start. \n", me);
      fprintf( stderr, " me = %ld clustrinv,mgs cannot handle this. \n", me);
      exit( -1 );
    }

    if( send_num == 0 ) {
      itype = 99;

      c1     = schedule[4*recv_cl];
      csiz   = schedule[4*recv_cl+1] - c1 + 1;
      
      blksiz = schedule[4*recv_cl+3] - schedule[4*recv_cl+2] + 1;

      nvecs  = count_list( recv_from, &mapZ[c1], &csiz);

      isize  = sizeof( DoublePrecision ) * blksiz * nvecs;

      first_buf = dscrat;
      dscrat += nvecs * blksiz;

#ifdef DEBUG1
  fprintf(stderr, " me = %d Just before mxread isize = %d nvecs = %d \n", me, isize, nvecs );
#endif
  
      ival = mxread_( first_buf, &isize, &recv_from, &itype );

#ifdef DEBUG1
  fprintf(stderr, " me = %d Just after mxread \n", me );
  for( j = 0; j < nvecs*blksiz; j++)
       fprintf(stderr, " me = %d first_buf[%d] = %g \n",
                     me, j, first_buf[j]);
#endif
  
    }
  }

  cl_ptr = schedule;
  for (clustr_ptr= 0;  clustr_ptr < cl_num ; clustr_ptr++) {

    if( clustr_ptr == 0 && send_num > 0 ) {

        c1  = schedule[ 4*send_cl ];
        cn  = schedule[ 4*send_cl + 1 ];
        bb1 = schedule[ 4*send_cl + 2 ];
        bn  = schedule[ 4*send_cl + 3 ];

        if( clustr_ptr == send_cl )
          cl_ptr += 4;
    }
    else {

      c1 = *(cl_ptr++);
      cn = *(cl_ptr++);
      bb1 = *(cl_ptr++);
      bn = *(cl_ptr++);

      if(  send_num > 0 ) {
	if( clustr_ptr == send_cl ) cl_ptr += 4;
      }
    }


    if ( c1 < *imin ) {
      Zvec = 0;
    }
    else {
      Zvec = c1 - *imin;
    }
    
    
    blksiz = bn - bb1 + 1;
    csiz = cn - c1 + 1;
    onenrm = fabs(d[bb1]) + fabs(e[bb1 + 1]);
    tmp = fabs(d[bn])+ fabs(e[bn]);
    onenrm = max(onenrm,tmp);
    
    for (i = bb1 + 1; i < bn; ++i){
      tmp = fabs(d[i])+fabs(e[i])+fabs(e[i + 1]);
      onenrm = max(onenrm, tmp );
    }
    
    stpcrt = 5.*sqrt((DoublePrecision ) 1.0e-1 / (DoublePrecision ) blksiz);
    
    indx = 0;
    if ( csiz > 1 ){
      for (j = c1; j <= cn; j++ ){
	if ( mapZ[j] == me ) {
	  i = indx + Zvec;
	  mapvecZ[ i ] = j;
	  
	  /*
	    Initialize vector to zero.
	    */
	  
	  for ( jjj = 0; jjj < bb1; jjj++ )
	    vecZ[i][jjj] = 0.0e0;
	  
	  /*
	    fill with random entries
	    */
	  
	  dlarnv_(&IONE, iseed, &blksiz, &vecZ[i][bb1]);
	  
          for ( jjj = bn; jjj < msize; jjj ++ )
	    vecZ[i][jjj] = 0.0e0;
	  
	  indx++;
	}
	else {
	  dlarnv_(&IONE, iseed, &blksiz, dscrat);
	}
      }
    }
    
#ifdef DEBUG1
    fprintf(stderr, " just before inv_it and mgs of clustrxx me = %d c1 = %d cn = %d \n", me, c1, cn );
#endif
    
    first = 0;
    if( clustr_ptr == 0 && send_num > 0 )
      first = 1;
    
    /*
      fine problem
      */
    
    if ( csiz == 1 ){
      msize = *n;
      for ( jjj = 0; jjj < bb1; jjj++ )
	vecZ[Zvec][jjj] = 0.0;
      blksiz = bn - bb1 + 1;
      dlarnv_(&IONE, iseed, &blksiz, &vecZ[Zvec][bb1]);
      /*
	for ( jjj = bb1; jjj < bn; jjj++ )
	vecZ[Zvec][jjj] = 1.0;
	*/
      for ( jjj = bn; jjj < msize; jjj++ )
	vecZ[Zvec][jjj] = 0.0;
      
      itmp = inv_it3( n, &c1, &cn, &bb1, &bn, &Zvec, mapZ, mapvecZ, vecZ, d, e, eval, &eps, &stpcrt, &onenrm, iscratch, dscrat);
      
      /*
	itmp = inv_it5( n, &c1, &cn, &bb1, &bn, &Zvec, mapZ, mapvecZ, vecZ,
	dplus, lplus, ld, lld, eval, &eps, &stpcrt,\
	&onenrm, iscratch, dscrat);
*/

    }
    else {
      /*
	fine cluster
      */
      
      itime = 2;
      
      for ( j = 0; j < INV_TIME; j++ ) {
	/*
	  itmp = inv_it4( &j, n, &c1, &cn, &bb1, &bn, &Zvec, mapZ, mapvecZ, vecZ, dplus, lplus, ld, lld, eval, &eps, &stpcrt, &onenrm, iscratch, dscrat);
	*/
	
	  itmp = inv_it3( n, &c1, &cn, &bb1, &bn, &Zvec, mapZ, mapvecZ, vecZ, d, e, eval, &eps, &stpcrt, &onenrm, iscratch, dscrat);
	
	if( itmp >  0 )
	  if( ibad == 0 || itmp < ibad ) 
	    ibad = itmp;
	
	for ( i = 0; i < itime ; i++ )
	  mgs_3( &csiz, vecZ, &mapZ[c1], &bb1, &bn, &Zvec, &first, first_buf, iscratch, dscrat);
	
	itime = itime-1 ;
      }
    }
    
#ifdef DEBUG1
    fprintf(stderr, " clustrxx3 me = %d before send/rec \n", me );
#endif
    
    /*
     * Swap beginning portions of clusters which are distributed
     * across more than one processor.
     */

  if( clustr_ptr == 0 && send_num > 0 ) {
    
      itype = 99;

      if( recv_num > 0  &&  ( myindx % 2 ) == 0 ) { 
        xc1     = schedule[4*recv_cl];
        xcsiz   = schedule[4*recv_cl+1] - xc1 + 1;

        xblksiz = schedule[4*recv_cl+3] - schedule[4*recv_cl+2] + 1;
	
        nvecs  = count_list( recv_from, &mapZ[xc1], &xcsiz);
        isize = sizeof( DoublePrecision ) * xblksiz * nvecs;
	
        first_buf = dscrat;
        dscrat += nvecs * xblksiz;

#ifdef DEBUG1
	fprintf(stderr, " me = %d Just before mxread 2 isize = %d nvecs = %d \n", me, isize, nvecs );
#endif
	
        ival = mxread_( first_buf, &isize, &recv_from, &itype );
      }
      
      nvecs = 0;
      for (j = c1; j <= cn; j++ ){
        if ( mapZ[j] == me ) {
	  dcopy_(&blksiz, &vecZ[nvecs+Zvec][bb1], &IONE,
                 dscrat+nvecs*blksiz, &IONE );
          nvecs++;
        }
      }
	
      isize = sizeof( DoublePrecision ) * blksiz * nvecs;
#ifdef DEBUG1
  fprintf(stderr, " me = %d Just before mxwrit isize = %d nvecs = %d \n", me, isize, nvecs );
  for( j = 0; j < nvecs*blksiz; j++)
       fprintf(stderr, " me = %d sending dscrat[%d] = %g \n",
                     me, j, dscrat[j]);
#endif
  
      ival = mxwrit_( dscrat, &isize, &send_to, &itype );

      if( recv_num > 0  &&  ( myindx % 2 ) != 0 ) {
        xc1     = schedule[4*recv_cl];
        xcsiz   = schedule[4*recv_cl+1] - xc1 + 1;
        xblksiz = schedule[4*recv_cl+3] - schedule[4*recv_cl+2] + 1;
	
        nvecs  = count_list( recv_from, &mapZ[xc1], &xcsiz);
        isize = sizeof( DoublePrecision ) * xblksiz * nvecs;
	
        first_buf = dscrat;
        dscrat += nvecs * xblksiz;

#ifdef DEBUG1
  fprintf(stderr, " me = %d Just before mxread 3 isize = %d nvecs = %d \n", me, isize, nvecs );
#endif
  
        ival = mxread_( first_buf, &isize, &recv_from, &itype );
      }
    }
#ifdef DEBUG1
  fprintf(stderr, " clustrxx3 me = %d after send/rec \n", me );
#endif

  }
  
#ifdef DEBUG1
  fprintf(stderr, " me = %d Exiting clustrinv_ \n", me );
#endif
  
  return(ibad);
}
