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
 *  -- PEIGS  routine (version 2.2) --
 *     Pacific Northwest Laboratory
 *     July 28, 1995
 *
 *======================================================================
 */
#include <stdio.h>
#include <memory.h>
#include <math.h>

#include "blas_lapack.h"
#include "globalp.c.h"

#define PANELSIZE 5

#define max(a,b) ((a) > (b) ? (a) : (b))
#define min(a,b) ((a) < (b) ? (a) : (b))
#define ffabs(a) ((a) > (0.) ? (a) : (-a))
#define sgn(a) ((a) >= (0.) ? (1.e0) : (-1.e0))

/*
   Internal PeIGS routine
   
   modified gram-schmidt, 1-D systolic wrap panelized version
   */

void mgs_3( n, colF, mapF, b1, bn, nvecsZ, first, first_buf, iscratch, scratch)
     Integer *n, mapF[], *b1, *bn, *nvecsZ, *first, iscratch[];
     DoublePrecision **colF, first_buf[], scratch[];
{
  
  /*
     n = number of vector to be orthogonalized
     b1 = beginning of vector subscript
     bn = end of vector subscript (e.g. colF[i][b1:bn] )
     colF = double pointer to the matrix
     mapF = array describing the processors holding columns of F
     first = 1, then mgs only my local block
                otherwise do full mgs.
     iscratch = integer scratch space
     scratch = double precision scratch space
     */
  
  static Integer IONE = 1, MONE = (DoublePrecision) -1.0e0;
  Integer jndx, vec_len;
  Integer i, k, me, isize, indx;
  Integer nvecs_in, nvecs, iii;
  Integer bb, nproc;
  Integer rsize, mvecs, kk, iter;
  
  Integer *nleft;
  Integer *iscrat, *proclist, naproc;
  Integer me_indx, max_vecs, max_panel;

  DoublePrecision *buffer, *in_buffer;
  DoublePrecision t, *dptr, *dptr1;

  /* blas calls */

  extern void dcopy_(), daxpy_(), dscal_();
  extern DoublePrecision ddot_();
  extern DoublePrecision gdot_();
  extern DoublePrecision dnrm2_();

  extern void bbcast00();
  extern Integer count_list();
  extern Integer mxmynd_();
  extern Integer reduce_list4();
  extern Integer indxL();
  
  
  me = mxmynd_();
  naproc = mxnprc_();
  
#ifdef DEBUG1
   fprintf(stderr, " me = %d in mgs \n", me );
#endif
#ifdef DEBUG2
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
     fprintf(stderr, " mgs5 me = %d mapZ[%d] = %d \n", me, iii, mapF[iii]);
#endif
   
   nvecs = count_list( me, mapF, n );
   
   if ( nvecs == 0 )
     return;
   
   vec_len = *bn - *b1 + 1;
   bb      = *b1;
   
   iscrat = iscratch;
   
   proclist = iscrat;
   nproc = reduce_list4( *n, mapF, proclist, &iscrat[naproc] );
   iscrat += nproc;
   
   nleft = iscrat;
   iscrat += nproc;
   
#ifdef DEBUG1
   fprintf(stderr, " me = %d mgs nprocs = %d nvecs = %d \n", me , nproc, nvecs);
#endif

   if( *first == 1  || nproc == 1 ) {
     
     /* mgs local block and return */
     
     k = *nvecsZ;
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
     return;
   }
   
   k = *nvecsZ;
   
   buffer = (DoublePrecision *) scratch;
   in_buffer = buffer;
   
   /*
    *  MGS my part of cluster against the first part of the cluster,
    *  which is stored in first_buf.
    */
   
   nvecs_in = count_list( proclist[0], mapF, n );
   
   if( me == proclist[1] ) {
     isize = vec_len * nvecs_in;
     dcopy_( &isize, first_buf, &IONE, in_buffer, &IONE);
   }
   
   if( nproc > 2 ) {
     rsize = nvecs_in * vec_len * sizeof(DoublePrecision);
     
     if ( rsize != 0 )
       bbcast00( (char *) in_buffer, rsize, 999,  proclist[1],
		 nproc-1, &proclist[1] );    
   }
   
   dptr = in_buffer;
   for ( iii = k; iii < k + nvecs; iii++ ){
     dptr = &colF[iii][bb];
     dptr1 = in_buffer;
     for ( jndx = 0; jndx < nvecs_in; jndx++ ){
       t = ddot_( &vec_len, dptr, &IONE, dptr1, &IONE);
       t *= MONE;
       daxpy_( &vec_len, &t, dptr1, &IONE, dptr, &IONE );
       dptr1 += vec_len;
     }
   }
   
   /* 
    * MGS the rest of the cluster.
    */
   
   me_indx = indxL( me, nproc, proclist);
   
   /* 
    * Determine maximum number of panels owned by
    * any processor in this cluster
    */
   
   max_vecs = 0;
   for ( i = 1; i < nproc; i++ ) {
     kk       = count_list( proclist[i], mapF, n );
     nleft[i] = kk;
     max_vecs = max( max_vecs, kk );
   }
   
   max_panel = max_vecs / PANELSIZE;
   if( max_panel * PANELSIZE != max_vecs )
     max_panel++;
   
   kk = k;
   for ( iter = 0; iter < max_panel; iter++ ){
     for( i = 1; i < nproc; i++ ) {
       mvecs = min( nleft[i], PANELSIZE );
       if ( mvecs == 0 ) 
	 continue;
       nleft[i] -= mvecs;
       if ( me_indx == i ) {
        for ( jndx = kk; jndx < kk + mvecs; jndx++ ){
 	  dptr = &colF[jndx][bb];
 	  t = dnrm2_( &vec_len, dptr, &IONE );
	if ( t != 1.0e0 ) {
	  t = (DoublePrecision) 1.0e0/t;
	  dscal_( &vec_len, &t, dptr, &IONE);
	}
	  for ( indx = jndx + 1; indx < kk + mvecs; indx++ ){
	    dptr1 = &colF[indx][bb];
	    t = ddot_( &vec_len, dptr, &IONE, dptr1, &IONE );
	   if ( fabs(t) > DLAMCHE ) {
	    t *= MONE;
	    daxpy_( &vec_len, &t, dptr, &IONE, dptr1, &IONE );
	}
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
      
      rsize = mvecs * vec_len * sizeof(DoublePrecision);
      
      bbcast00( (char *) in_buffer, rsize, 888,  proclist[i],
                nproc-1, &proclist[1] );    
      
      /*
       * the buffer contains incoming orthonormal vectors
       */
      
      dptr = in_buffer;
      for ( iii = kk; iii < k + nvecs; iii++ ){
        dptr = &colF[iii][bb];
        dptr1 = in_buffer;
        for ( jndx = 0; jndx < mvecs; jndx++ ){
          t = ddot_( &vec_len, dptr, &IONE, dptr1, &IONE);
	if ( fabs(t) > DLAMCHE ) {
	  t *= MONE;
	  daxpy_( &vec_len, &t, dptr1, &IONE, dptr, &IONE );
	}
          dptr1 += vec_len;
	}
      }
     }
   }
   return;
}


