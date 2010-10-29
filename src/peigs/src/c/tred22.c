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
/************************************************
 *
 *     C routine for reduction of a symmetric
 *     matrix to symmetric tridiagonal form
 *
 *  the tridiagonal matrix is returned as d[0:n-1]
 *  and e[1:n-1]
 *
 * Integer tred2 (n, vecA, mapA, Q, mapQ, diag, upperdiag, iwork, work )
 *     DoublePrecision **vecA,               matrix to be reduced
 *       **Q,                       eigenvector matrix
 *       *diag,              diagonal elements of tri-diagonal matrix
 *       *upperdiag,         upper diagonal elements of tri-diagonal matrix
 *       *work;              Householder vector (temp space of size n doubles)
 *     
 *       Integer   *n,                  problem size
 *       *mapA,
 *       *mapQ,
 *       *iwork;          integer scratch work space
 *
 *
 *
 *
 */

/* ================================================================= */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "globalp.c.h"

extern DoublePrecision ddot_(), dnrm2_(), dasum_();

#define MSG_START 25000

/*
  ----- interfaces to Littlefield's portable mx... comm routines -------- 
  */

/* 
  The following #define's are to avoid recursive calls if the mx...
  routines are in turn implemented using PICL.
  */

#define sync0  prs1sync0
#define clock0 prs1clock0
#define who0   prs1who0

#define MAXPROCS 1024
extern DoublePrecision clock0();

void bbcast00(buf, len, type, root, snumprocs, plist)
     char *buf;
     Integer len;
     Integer type;
     Integer root, snumprocs;
     Integer *plist;
{
  /*
    mkplist ();  assume that plist has the list of all the processors
    */

  extern Integer mxbrod_();
  
  mxbrod_(buf,&root,&len,&snumprocs,plist,&type);
}

void ibcast00(buf, len, type, root, snumprocs, plist)
     Integer *buf;
     Integer len;
     Integer type;
     Integer root, snumprocs;
     Integer *plist;
{
  /*
    mkplist ();  assume that plist has the list of all the processors
    */
  
  extern Integer mxbrod_();
  
  mxbrod_(buf,&root,&len,&snumprocs,plist,&type);
}

void peigs_gmax00(buf, items, datatype, msgtype, root, snumprocs, plist, work)
     char *buf;
     Integer items;
     Integer datatype;
     Integer msgtype;
     Integer root;
     Integer *plist, snumprocs;
     DoublePrecision *work;  /* workspace containing at least bufsiz bytes (see cmbbrf.h) */
{
  Integer eight = sizeof(DoublePrecision);
  extern Integer maxd_();
  extern Integer mxcombv1_();
  
  mxcombv1_( buf, maxd_, &eight, &items, &snumprocs, plist, &msgtype, (char *)work);
}




void gsum00(buf, items, datatype, msgtype, root, snumprocs, plist, work)
     /*
	Note that this implementation of gsum00 differs from Chinchalkar's
	in that it leaves the result in all nodes participating nodes
	(0..snumprocs-1) instead of just the root.  This allows deleting
	the explicit bbcast00 calls following each of the gsum00's.
	*/
     char *buf;
     Integer items;
     Integer datatype;
     Integer msgtype;
     Integer root;
     Integer *plist, snumprocs;
     DoublePrecision *work;  /* workspace containing at least bufsiz bytes (see cmbbrf.h) */
{
  Integer eight = sizeof(DoublePrecision);
  extern Integer sumdv_();
  extern Integer mxcombv1_();
  
  mxcombv1_( buf, sumdv_, &eight, &items, &snumprocs, plist, &msgtype, (char *)work);
}

void gsum01(buf, items, datatype, msgtype, root, snumprocs, plist, work)
     /*
       vectorized integer global sum
       */
     char *buf;
     Integer items;
     Integer datatype;
     Integer msgtype;
     Integer root;
     Integer *plist, snumprocs;
     DoublePrecision *work;  /* workspace containing at least bufsiz bytes (see cmbbrf.h) */
{
  Integer data_len = sizeof(Integer);
  extern Integer sumiv_();
  extern Integer mxcombv1_();
  
  mxcombv1_( buf, sumiv_, &data_len, &items, &snumprocs, plist, &msgtype, (char *)work);
}

/* ----- FORTRAN interface ---------- */


Integer tred2(n, vecA, mapA, Q, mapQ, diag, upperdiag, iwork, work )
     Integer   *n,                  /* problem size */
     *mapA,
     *mapQ,
     *iwork;          /* integer scratch work space */
     
     DoublePrecision **vecA,               /* matrix to be reduced */
     **Q,                /* eigenvector matrix */
     *diag,             /* diagonal elements of tri-diagonal matrix */
     *upperdiag,        /* upper diagonal elements of tri-diagonal matrix */
     *work;           /* Householder vector plus gsum00 workspace.
                         (temp space of size n +1 doubles) PLUS
                         at least bufsiz bytes (see cmbbrf.h) */
{
  
  /*
   * TRED2 : reduction of a real symmetric matrix to tri-diagonal form
   *         Uses Householder reductions.
   *
   * INPUT:
   *
   * A[n][n]         ->  input matrix to be reduced.
   *                     only the upper triangular portion need be filled
   *                     distributed by rows(or columns) in a wrap fashion. (not compact)
   * n               ->  problem size.
   * numprocs        ->  number of processors to use for solution
   *
   * OUTPUT:
   *
   * diag[n]         ->  diagonal elements of the tri-diagonal matrix
   *                     obtained from A
   * upperdiag[n]    ->  upper diagonal elements of the tri-diagonal matrix
   *                     obtained from A. upperdiag[0] is set to 0.0 and
   *                     upperdiag[1] to upperdiag[n-1] contains the desired
   *                     result
   * Q(n,n)          ->  eigenvector matrix
   *                     distributed by rows(or columns) in a wrap fashion (compact)
   *                      P_0.P_1.P_2...
   *
   * tred2() returns an error code, which, if 0 indicates success. NYI
   
   March 31, 1993
   modified to work with mapping list wrapping
   mapQ is assume to be the same as mapA
   
   */
  
  static Integer IONE = 1;
  DoublePrecision DZERO = (DoublePrecision) 0.0e0, DONE = (DoublePrecision) 1.0e0;
  
  Integer i, j, k;               /* counters */
  Integer row_indx, linfo;
  Integer msize;
  Integer n_procs, *mapvecA, *mapvecQ;
  Integer *iscrat, *proclist;

  DoublePrecision norm2,              /* 2 norm of a vector */
  new_norm,           /* 2 norm of a vector */
  first,              /* first component of a vector */
  onenorm,            /* one norm of a vector */
  p_t_v;              /* p'v */
  
  DoublePrecision w, temp_factor;     /* temp vars */
  DoublePrecision *HH_vec, *workMX;

  char msg[45];
  
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

  
  Integer size, nrowsA, nrowsQ;
  
  Integer iii;

  /* ------------------------------------------------------------------- */

  me = mxmynd_();
  strcpy ( msg, "Error in TRED2 \n");
  
  linfo = 0;
  
  if ( n == NULL ) {
    linfo = -1;
    xerbla_( "TRED2 \n", &linfo);
    return(1);
  }
  
  if ( *n < 1) {
    linfo = -1;
    xerbla_( "TRED2 \n", &linfo);
    return(1);
  }
  
  if ( vecA == NULL ) {
    linfo = -2;
    xerbla_( "TRED2 \n", &linfo);
    return(1);
  }
  
  if ( mapA == NULL ) {
    linfo = -3;
    xerbla_( "TRED2 \n", &linfo);
    return(1);
  }
  
  
  iii = mxnprc_();
  iscrat = mapA;
  for ( i = 0; i < *n; i++ ) {
    j = *(iscrat++);
    if ( j < 0 ) {
      linfo = -3;
      xerbla_( "TRED2 \n", &linfo);
    }
    if ( j > iii  ){
      linfo = -3;
      xerbla_( "TRED2 \n", &linfo);
    }
  }

  nrowsA = count_list( me, mapA, n);
  for ( i = 0 ; i < nrowsA; i++ )
    if ( vecA[i] == NULL ){
      linfo = -2; 
      fprintf(stderr, "node = %d NULL vector assignment in vecA \n", me);
      xerbla_( "TRED2 \n", &linfo);
    }
    
  if ( iwork == NULL ) {
    linfo = -8;
    xerbla_( "TRED2 \n", &linfo);
  }
  
  if ( work == NULL ){
    linfo = -9;
    xerbla_( "TRED2 \n", &linfo);
  }
  
  if ( Q == NULL )
    linfo = -4;
  
  if ( mapQ == NULL )
    linfo = -5;
  
  if ( diag == NULL )
    linfo = -6;

  if ( upperdiag == NULL )
    linfo = -7;
  
  sprintf(msg, "TRED22:Error in argument %d \n", linfo);

  /*
  printf(" in tred22.c me = %d \n", me );
  */
  
  
  g_exit_( &linfo, msg, mapA, n, iwork, work);
  
  /* checking the other arguments to see if further error check is possible */
  
  *iwork = *n;
  iscrat = iwork + 1;
  for ( k = 0; k < *n; k++ )
    *(iscrat++) = mapA[k];
  
  linfo = 0;
  k = (*n + 1)*sizeof(Integer);
  iscrat += *n + 1;
  pxerbla2_( &k, (char *)iwork, mapA, n, iscrat, &linfo );
  if( linfo != 0 )
    linfo = -1;

  g_exit_( &linfo, "TRED22: Mapping inconsistancies. mapA\n", mapA, n, iwork, work);
  
  iscrat = iwork;
  linfo = 0;
  k = *n * sizeof(Integer);
  pxerbla2_( &k, (char *) mapQ, mapA, n, iscrat, &linfo );
  
  mapdif1_( n, mapQ, n, mapA, iscrat, &j );
  if( linfo != 0  || j != 0 )
    j = -1;
  
  g_exit_( &j, "TRED22: Mapping set differ:mapQ and mapA\n", mapA, n, iwork, work);
  
  
  
  /*
    setting the location of the Q matrix
    */
  

  
  /*
    
    set the matrix Q to be distributed as matrix A  ; replicate this information
    One should probably do a more general Q matrix distribution; this is just to enforce
    that mapQ = mapA
    
    */
  
  /*
    DZERO out upperdiag[0]
    */
  
  upperdiag[0] = DZERO;
  
  /*
    set Q to the identity matrix
    */
  
  if ( n == NULL )
    exit(-1);
  
  
  if ( *n < 0 )
    exit(-1);
  
  msize = *n;
  
  HH_vec = work;

  workMX = work + *n + 1;

  iscrat = iwork;
  mapvecQ = iscrat;
  nrowsQ = fil_mapvec_( &me, &msize, mapQ, mapvecQ );
  iscrat += nrowsQ;
  mapvecA = iscrat;
  nrowsA = fil_mapvec_( &me, &msize, mapA, mapvecA );
  iscrat += nrowsA;
  proclist = iscrat;
  n_procs = reduce_list2( msize, mapA, proclist);  
  iscrat += n_procs;
  
  /*
    initialize Q matrix to the identity matrix
    */
  
  k = 0;
  for (i=0; i< msize; i++) {
    if ( mapQ[i] == me ) {
      for (j=0; j< msize; j++) {
	Q[k][j] = DZERO;
      }
      Q[k][i] = DONE;
      k++;
    }
  }
  
  
  /* apply the n-2 HH transformations */

  row_indx = -1;
  for (i=0; i<(msize-2); i++) {
    if ( mapA[i] == me ) {
      row_indx++;
      /*
	construct the HH vector based on ROW i
	*/
      
      /*
	scale the vector if necessary
	*/
      
      size = msize - i - 1;
      onenorm = dasum_(&size, &vecA[row_indx][1], &IONE);
      
      /* now we know the diagonal entry */
      
      diag[i] = vecA[row_indx][0];
      
      if (onenorm == DZERO) {
	/* then skip this transformation */
	
	vecA[row_indx][0] = DZERO;
	upperdiag[i+1] = DZERO;
      }
      else {
	vecA[row_indx][0] = DONE;
	
	/* scale vector */
	
	size = msize - i - 1;
	norm2 = dnrm2_(&size, &vecA[row_indx][1], &IONE)/onenorm;
	vecA[row_indx][1] /= onenorm;
	
	if (vecA[row_indx][1] < DZERO)  {
	  first = vecA[row_indx][1] - norm2;
	  upperdiag[i+1] = norm2*onenorm;
	}
	else {
	  first = vecA[row_indx][1] + norm2;
	  upperdiag[i+1] = -norm2*onenorm;
	}
	
	/* change A[i][i+1:n-1] to v and normalize it too */
	
	new_norm = norm2*norm2 - vecA[row_indx][1] * vecA[row_indx][1] + first*first;
	new_norm = sqrt(new_norm);
	
	vecA[row_indx][1] = first/new_norm;
	temp_factor = DONE/onenorm/new_norm;
	size = msize - i - 2;
	dscal_(&size, &temp_factor, &vecA[row_indx][2], &IONE);
      }
      
      /* send v to all other nodes */
      
      bbcast00( (char * ) &vecA[row_indx][0], (msize-i)*sizeof(DoublePrecision), i, me, n_procs, proclist);
      
      /* copy the new HH vector into HH_vec */
      
      size = msize - i;
      dcopy_(&size, &vecA[row_indx][0], &IONE, &HH_vec[i], &IONE);
      
    }
    else  {
      /*
	ireceive v
	*/
      bbcast00( (char *) &HH_vec[i], (msize-i)*sizeof(DoublePrecision), i, mapA[i], n_procs, proclist);
      diag[i] = DZERO;
      upperdiag[i+1] = DZERO;
    }
    
    if (HH_vec[i] != DZERO ) {
      /* then this transformation has not been skipped */
      
      /*
       * Now HH_vec[i+1:n-1] = v and v'v = 1
       * update A(i+1:n-1, i+1:n-1)
       * NOTE: diag[i+1:n-1] is unused. Store p in it.
       */
      
      /* DZERO out vector; should dscal it */
      
      for (k=i+1; k<msize; k++) {
	diag[k] = DZERO;
      }
      
      j = row_indx;
      for ( k = i + 1; k < msize; k++) {
	if ( mapA[k] == me ) {
	  size = msize - k - 1;
	  j++;
	  
	  /* above the diagonal */
	  
	  diag[k] += ddot_( &size, &vecA[j][1], &IONE, &HH_vec[k+1], &IONE);
	  
	  /* below the diagonal */
	  
	  daxpy_(&size, &HH_vec[k], &vecA[j][1], &IONE, &diag[k+1], &IONE);
	  
	  /* diagonal entry */
	  
	  diag[k] += vecA[j][0] * HH_vec[k];
	}
      }
      
      
      /* p is now distributed across processors. collect it */
      
      gsum00( (char *) &diag[i+1], msize-i-1, 5, MSG_START, mapA[0], n_procs, proclist, workMX);
      
      size = msize - i - 1;
      temp_factor = ((DoublePrecision) -2.0e0);
      dscal_(&size, &temp_factor, &diag[i+1], &IONE );
      
      /* over-write diag[i+1:n-1] with w */
      
      p_t_v = -ddot_(&size, &diag[i+1], &IONE, &HH_vec[i+1], &IONE );
      daxpy_(&size, &p_t_v, &HH_vec[i+1], &IONE, &diag[i+1], &IONE);
      
      /* and finally update A[i+1:n-1][i+1:n-1] */
      
      j = row_indx;
      for (k=i+1; k<msize; k++) {
	if ( mapA[k] == me) {
	  j++;
	  size = msize - k;
	  daxpy_(&size, &HH_vec[k], &diag[k], &IONE , &vecA[j][0], &IONE );
	  daxpy_(&size, &diag[k], &HH_vec[k], &IONE , &vecA[j][0], &IONE );
	}
      }
      
      /* FORWARD ACCUMULATION */
      /* compute w. No need to store it. beta is -2 */
      
      for (j=0; j < nrowsQ; j++) {
	size = msize - i - 1;
	w =((DoublePrecision) -2.0e0) * ddot_(&size, Q[j] + i + 1, &IONE , &HH_vec[i+1], &IONE );
	daxpy_(&size, &w, &HH_vec[i+1], &IONE, Q[j] + i + 1, &IONE);
      }
    }
  }
  
  /*
    now whatever is left over
    */
  
  if (msize == 1) {
    if (mapA[0] == me ) {
      diag[0] = vecA[0][0];
    }
  }
  else {
    if ( mapA[msize-2] == me  ){
      j = indxL( msize-2, nrowsA, mapvecA);
      upperdiag[msize-1] = vecA[j][1];
      diag[msize-2] = vecA[j][0];
    }
    else {
      upperdiag[msize-1] = DZERO;
      diag[msize-2] = DZERO;
    }
    
    if ( mapA[msize-1] == me ){
      j = indxL( msize-1, nrowsA, mapvecA);
      diag[msize-1] = vecA[j][0];
    }
    else {
      diag[msize-1] = DZERO;
    }
  }
  
  /*
    collect diag[] and upperdiag[]
    */
  
  gsum00((char * ) diag, msize, 5, MSG_START+2, mapA[0], n_procs, proclist, workMX);
  gsum00((char * ) &upperdiag[1], msize-1, 5, MSG_START+4, mapA[0], n_procs, proclist, workMX);

  
  /*
    printf(" out tred22.c me = %d \n", me );
    */

  /*
  for (k=0; k<msize; k++)
    printf(" k = %d tridiag d = %f e = %f \n", k, diag[k], upperdiag[k]);
    */

  return 0;
}
