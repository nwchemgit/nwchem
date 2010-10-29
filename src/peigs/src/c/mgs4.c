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
#include <memory.h>
#include <math.h>

#include "globalp.c.h"

#define PANELSIZE 10

#define max(a,b) ((a) > (b) ? (a) : (b))
#define min(a,b) ((a) < (b) ? (a) : (b))

/*
   Internal PeIGS routine
   
   modified gram-schmidt, 1-D systolic panelized version
   */

void mgs_3( n, colF, mapF, b1, bn, nvecsZ, first, first_buf, iscratch, scratch)
     Integer *n, mapF[], *b1, *bn, *nvecsZ, *first, iscratch[];
     DoublePrecision **colF, first_buf[], scratch[];
{
  /*
   */
  
  /*
     n = number of vector to be orthogonalized
     b1 = beginning of vector subscript
     bn = end of vector subscript
     (e.g. colF[i][b1:bn] )
     colF = double pointer to the matrix
     mapF = array describing the processors holding columns of F
     first = 1, then mgs only my local block
                otherwise do full mgs.
     iscratch = integer scratch space
     scratch = double precision scratch space
     */
  
  static Integer IONE = 1, MONE=(DoublePrecision) -1.0e0;
  Integer jndx, vec_len;
  Integer i, k, me, isize, indx, kndx;
  Integer nvecs_in, nvecs, iii;
  Integer j, bb, nproc;
  Integer rsize, miter, mvecs, kk, itype, iter, iremain;
  
  Integer *mapvecF;
  Integer *iscrat, *proclist, naproc;
  Integer *mapvec_in, me_indx;

  DoublePrecision *buffer, *in_buffer;
  DoublePrecision t, *dptr, *dptr1;

  /*
    blas calls
    */

  extern void dcopy_(), daxpy_(), dscal_();
  extern DoublePrecision ddot_();
  extern DoublePrecision dnrm2_();
  extern void bbcast00();

  /* 
    mxsubs calls
    */

  extern Integer mxwrit_(),  mxread_();
  extern Integer count_list();
  extern Integer fil_mapvec_();
  extern Integer mxmynd_();
  extern Integer reduce_list4();
  extern Integer indxL();
  
  /*
    at this point inputs are minimally acceptable

    check to see if mapF are the same set of processors
    */

  me = mxmynd_();
  naproc = mxnprc_();


#ifdef DEBUG1
  fprintf(stderr, " \n" );
  fprintf(stderr, " In mgs1b me = %d \n", me );
  fprintf(stderr, " \n" );
  i = *nvecsZ-1;
  for( iii = 0; iii < *n; iii++)
    if( mapF[iii] == me ) {
      i++;
      for( j = *b1; j <= *bn; j++)
       fprintf(stderr, " mgs1b me = %d vecZ[%d][%d] = %g \n",
                     me, iii, j, colF[i][j]);
    }

  fprintf(stderr, " \n" );
  for( iii = 0; iii < *n; iii++)
  fprintf(stderr, " mgs1b me = %d mapZ[%d] = %d \n", me, iii, mapF[iii]);
  fprintf(stderr, " in mgs1b me = %d \n", me);
#endif

  vec_len = *bn - *b1 + 1;

  iscrat = iscratch;

  mapvecF = iscrat;
  nvecs = fil_mapvec_( &me, n, mapF, mapvecF );
  iscrat += nvecs;
  
  if ( nvecs == 0 )
    return;
  
  bb = *b1;
  
  proclist = iscrat;
  nproc = reduce_list4( *n, mapF, proclist, iscrat + naproc );

#ifdef DEBUG1
  fprintf(stderr, " me = %d mgs1b nprocs = %d n = %d \n", me , nproc, *n);
#endif

  if( *first == 1  || nproc == 1 ) {
  
    /*
     * mgs local block and return
     */
      
    k = *nvecsZ;
    for ( jndx = k; jndx < k + nvecs; jndx++ ){
	dptr = &colF[jndx][bb];
	t = dnrm2_( &vec_len, dptr, &IONE );
	t = (DoublePrecision) 1.0e0/t;
	dscal_( &vec_len, &t, dptr, &IONE);
	for ( indx = jndx + 1; indx < k + nvecs; indx++ ){
	  dptr1 = &colF[indx][bb];
	  t = -ddot_( &vec_len, dptr, &IONE, dptr1, &IONE );
	  daxpy_( &vec_len, &t, dptr, &IONE, dptr1, &IONE );
	}
    }
    return;
  }
      
  iscrat += nproc;
  mapvec_in = iscrat;
  
  me_indx = indxL( me, nproc, proclist);
  
  buffer = (DoublePrecision *) scratch;
  in_buffer = buffer;

  k = *nvecsZ;
  
  /*
   *  MGS my part of cluster againt the first part of the cluster,
   *  which is stored in first_buf.
   */

  nvecs_in = fil_mapvec_( &proclist[0], n, mapF, mapvec_in);

  if( me == proclist[1] ) {
    isize = vec_len * nvecs_in;
    dcopy_( &isize, first_buf, &IONE, in_buffer, &IONE);
  }

  if( nproc > 2 ) {
    rsize = nvecs_in * vec_len * sizeof(DoublePrecision);
      
    if ( rsize != 0 )
      bbcast00( (char *) in_buffer, rsize, 11113,  proclist[1],
                nproc-1, &proclist[1] );    
  }

  dptr = in_buffer;
  for ( iii = k; iii < k + nvecs; iii++ ){
    dptr = &colF[iii][bb];
    dptr1 = in_buffer;
    for ( jndx = 0; jndx < nvecs_in; jndx++ ){
      t = -ddot_( &vec_len, dptr, &IONE, dptr1, &IONE);
      daxpy_( &vec_len, &t, dptr1, &IONE, dptr, &IONE );
      dptr1 += vec_len;
    }
  }

  /*
   *  MGS the rest of the cluster.
   */

  for ( i = 1; i < nproc - 1; i++ ) {
    
    if ( me_indx < i ) 
      break;
    
    kndx = proclist[i];

    nvecs_in = fil_mapvec_( &kndx, n, mapF, mapvec_in);

    miter   = nvecs_in / PANELSIZE;
    iremain = nvecs_in - miter * PANELSIZE;

    if( miter * PANELSIZE != nvecs_in )
       miter++;

    mvecs = PANELSIZE;
    rsize = mvecs * vec_len * sizeof(DoublePrecision);

    kk    = k;
    itype = 11113;
    for ( iter = 1; iter < miter + 1; iter++ ){

      if( iter == miter && iremain > 0 ) {
        mvecs = iremain;
        rsize = mvecs * vec_len * sizeof(DoublePrecision);
      }

      if ( kndx == me ) {
      
        for ( jndx = kk; jndx < kk + mvecs; jndx++ ){
 	  dptr = &colF[jndx][bb];
 	  t = dnrm2_( &vec_len, dptr, &IONE );
	  t = 1.0e0/t;
	  dscal_( &vec_len, &t, dptr, &IONE);
	  for ( indx = jndx + 1; indx < kk + mvecs; indx++ ){
	    dptr1 = &colF[indx][bb];
	    t = ddot_( &vec_len, dptr, &IONE, dptr1, &IONE );
	    t *= MONE;
	    daxpy_( &vec_len, &t, dptr, &IONE, dptr1, &IONE );
	  }
        }

        /* load up buffer */
      
   	dptr = in_buffer;
	for ( indx = kk; indx < kk + mvecs; indx++ ){
	  dcopy_( &vec_len, &colF[indx][bb], &IONE, dptr, &IONE);
	  dptr += vec_len;
        }

        kk += mvecs;
      }

      bbcast00( (char *) in_buffer, rsize, 11113,  proclist[i],
                nproc-i, &proclist[i] );    
      
      /*
       * the buffer contains orthonormal vectors
       */

      if ( iter < miter  || kndx != me ) {
        dptr = in_buffer;
        for ( iii = kk; iii < k + nvecs; iii++ ){
          dptr = &colF[iii][bb];
          dptr1 = in_buffer;
          for ( jndx = 0; jndx < mvecs; jndx++ ){
            t = ddot_( &vec_len, dptr, &IONE, dptr1, &IONE);
	    t *= MONE;
            daxpy_( &vec_len, &t, dptr1, &IONE, dptr, &IONE );
            dptr1 += vec_len;
          }
        }
      }
    }
  }

  /*
   * Last local mgs and no one to send to.
   */
      
  if( me_indx == nproc - 1 ) {

    for ( jndx = k; jndx < k + nvecs; jndx++ ){
      dptr = &colF[jndx][bb];
      t = dnrm2_( &vec_len, dptr, &IONE );
      t = 1.0e0/t;
      dscal_( &vec_len, &t, dptr, &IONE);
      for ( indx = jndx + 1; indx < k + nvecs; indx++ ){
        dptr1 = &colF[indx][bb];
        t = ddot_( &vec_len, dptr, &IONE, dptr1, &IONE );
	t *= MONE;
        daxpy_( &vec_len, &t, dptr, &IONE, dptr1, &IONE );
      }
    }
  }

#ifdef DEBUG2
  fprintf(stderr, " out mgs1b just before bbcast0 me = %d \n", me);
#endif
  
  /*
    making sure everyone in the local ring is synchronized
    */
  
 if ( nproc > 2 ) 
     bbcast00( (char *) &t, 1, 11112, proclist[nproc-1], nproc-1, &proclist[1]);

#ifdef DEBUG2
	fprintf(stderr, " really out of mgs1b me = %d \n", me );
  fprintf(stderr, " \n" );
#endif
#ifdef DEBUG1
  i = *nvecsZ-1;
  for( iii = 0; iii < *n; iii++)
    if( mapF[iii] == me ) {
      i++;
      for( j = *b1; j <= *bn; j++)
       fprintf(stderr, " mgs1b me = %d vecZ[%d][%d] = %g \n",
                     me, iii, j, colF[i][j]);

    }
  fprintf(stderr, " \n" );
  fprintf(stderr, " Out of mgs1b me = %d \n", me );
  fprintf(stderr, " \n" );
#endif
  
  return;
}
