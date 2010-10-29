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
#include <stdlib.h>

#include "globalp.c.h"

extern DoublePrecision ddot_(), dnrm2_(), dasum_();
extern void daxpy_(), dscal_();

#define MSG_START 25000

/*
  ----- interfaces to Littlefield's portable mx... comm routines -------- 
  */

/* 
  The following #define's are to avoid recursive calls if the mx...
  routines are in turn implemented using PICL.
  */

#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define SGN(a) ((a) > (0) ? (DoublePrecision) (1.0e0) : (DoublePrecision) (-1.0e0))
#define MIN(a,b) ((a) > (b) ? (b) : (a))
#define FABS(a) ((a) > ((DoublePrecision) 0.0e0 ) ? (a) : (-a))


#define sync0  prs1sync0
#define clock0 prs1clock0
#define who0   prs1who0

extern DoublePrecision clock0();

/* ----- FORTRAN interface ---------- */

Integer mgscs(n, vecA, mapA, b1, bn, c1, cn, iwork, work )
     Integer   *n,                  /* problem size */
     *mapA,
     c1, cn, b1, bn,
     *iwork;          /* integer scratch work space */
     
     DoublePrecision
     **vecA,               /* matrix to be reduced */
       *work;           /* Householder vector plus gsum00 workspace.
			   (temp space of size n +1 doubles) PLUS
                         at least bufsiz bytes (see cmbbrf.h) */
{
  
  /*
   * MGS : reduction of a real symmetric matrix to tri-diagonal form
   *         Uses Householder reductions.
   
   March 31, 1993
   modified to work with mapping list wrapping
   
   */
  
  static Integer IONE = 1;
  
  Integer i, j, k, ii;               /* counters */
  Integer column_indx, linfo;
  Integer msize;
  Integer n_procs, *mapvecA;
  Integer *iscrat, *proclist;

  DoublePrecision   onenorm;
  
  DoublePrecision t, *ptr, syncco[1];

  Integer me;                    /* my node number */
  
  
  char *strcpy();
  
  extern Integer mxmynd_();
  extern Integer mxmync_();
  extern Integer reduce_list2();
  extern Integer count_list();
  extern Integer indxL();
  extern Integer mxnprc_();
  extern Integer fil_mapvec_();

  extern void xerbla_();
  extern void mapdif1_();

  extern void g_exit_();
  extern void pxerbla2_();

  extern void dscal_();
  extern void dcopy_();
  extern void daxpy_();
  extern DoublePrecision dasum_();
  extern void bbcast00();

  
  Integer csize, ncolumnsA;
  
  Integer iii;

  /* ------------------------------------------------------------------- */

  me = mxmynd_();
  
  linfo = 0;
  
  if ( n == NULL ) {
    linfo = -1;
    xerbla_( "MGSCS \n", &linfo);
    return(1);
  }
  
  if ( *n < 1) {
    linfo = -1;
    xerbla_( "MGSCS \n", &linfo);
    return(1);
  }
  
  if ( vecA == NULL ) {
    linfo = -2;
    xerbla_( "MGSCS \n", &linfo);
    return(1);
  }
  
  if ( mapA == NULL ) {
    linfo = -3;
    xerbla_( "MGSCS \n", &linfo);
    return(1);
  }
  
  
  iii = mxnprc_();
  iscrat = mapA;
  for ( i = 0; i < *n; i++ ) {
    j = *(iscrat++);
    if ( j < 0 ) {
      linfo = -3;
      xerbla_( "MGSCS \n", &linfo);
    }
    if ( j > iii  ){
      linfo = -3;
      xerbla_( "MGSCS \n", &linfo);
    }
  }

  ncolumnsA = count_list( me, mapA, n);

  for ( i = 0 ; i < ncolumnsA; i++ )
    if ( vecA[i] == NULL ){
      linfo = -2; 
      fprintf(stderr, "node = %d NULL vector assignment in vecA \n", me);
      xerbla_( "MGSCS \n", &linfo);
    }
    
  if ( iwork == NULL ) {
    linfo = -8;
    xerbla_( "MGSCS \n", &linfo);
  }
  
  if ( work == NULL ){
    linfo = -9;
    xerbla_( "MGSCS \n", &linfo);
  }
  
  if ( n == NULL ) 
    exit(-1);
  
  if ( *n < 0 )
    exit(-1);
  
  msize = *n;
  ii = cn -c1 +1 ;
  iscrat = iwork;
  mapvecA = iscrat;
  ncolumnsA = fil_mapvec_( &me, &msize, mapA, mapvecA );
  iscrat += ncolumnsA;
  proclist = iscrat;
  n_procs = reduce_list2( ii, &mapA[c1], proclist);  
  iscrat += n_procs;
  
  ii = -1;
  for ( i = 0; i < cn+1; i++ ) {
    if ( mapA[i] == me ) {
      ii++;
      if ( i >= c1 )
	break;
    }
  }
  
  /*
     printf(" before me %d c1 %d cn %d b1 %d bn %d nprocs %d \n", me, c1, cn, b1, bn, n_procs);
     */
  
  if ( ii == -1 )
    return 0;
  
  /*  printf(" after me %d c1 %d cn %d b1 %d bn %d nprocs %d \n", me, c1, cn, b1, bn, n_procs);
   */
  
  
#ifdef DEBUG7
  printf(" me %d c1 %d cn %d b1 %d bn %d nprocs %d \n", me, c1, cn, b1, bn, n_procs);
  fflush(stdout);
#endif
  
  /*
     printf("n = %d  me %d c1 %d cn %d b1 %d bn %d nprocs %d \n", *n, me, c1, cn, b1, bn, n_procs);
     fflush(stdout);
     */
  
  
  column_indx = ii;
  csize = bn - b1 + 1 ;
  for (i=c1; i<= cn; i++) {
    if ( mapA[i] == me ) {
      ptr = &vecA[column_indx][b1];
      onenorm = dnrm2_(&csize, ptr, &IONE);
      t = 1.0e0/onenorm;
      dscal_( &csize, &t, ptr, &IONE);
      dcopy_( &csize, &vecA[column_indx][b1], &IONE, work, &IONE);
      column_indx++;
    }
    
    
    bbcast00( (char * ) work, (csize)*sizeof(DoublePrecision), c1,
	      mapA[i], n_procs, proclist);
    
    k = column_indx ;
    for ( j = i+1; j <= cn; j++ ){
      if ( mapA[j] == me ) {
	t = ddot_( &csize, work, &IONE, &vecA[k][b1], &IONE);
	t *= -1.0e0;
	daxpy_( &csize, &t, work, &IONE, &vecA[k][b1], &IONE);
	k++;
      }
    }
  }
  
  if ( n_procs > 1 ) {
    syncco[0] = 0.0e0;
    gsum00( (char *) syncco, 1, 0, 1, mapA[c1], n_procs, proclist, work);
  }
  
  return 0;
}


