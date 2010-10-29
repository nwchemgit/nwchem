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
#include <string.h>
#include <math.h>
#include <memory.h>
#include "globalp.c.h"

#define max(a,b) ((a) > (b) ? (a) : (b))
#define min(a,b) ((a) < (b) ? (a) : (b))

void pdcomplex(n, vecZ, mapZ, eval, scratch, iscratch, info)
     Integer            *n, *mapZ, *iscratch, *info;
     DoublePrecision    **vecZ, *eval, *scratch;
{
  
/*
 *
 *  Our parallel version of LAPACK's dspev.
 *
 *
 *  Purpose
 *  =======
 *
 reduce the complex hermitian 2*n real u + iv representation produced
from the peigs output to n real u + iv that are linear independent over i
 *
 *---------------------------------------------------------------------
 */
  
  /*
   * Local variables
   * ---------------
   */
  Integer             k, me, nn_proc, msize, irange, ivector, meigval, ilb,
                      iub, nproc, isize, nvecsA, nvecsZ, linfo, maxinfo,
                      i, j, itmp;
  Integer             *i_scrat, *proclist, lll, jjj;
  Integer l_good, nvecZ, cvecZ, IONE=1, ccvecZ, m;
  DoublePrecision     lb, ub, abstol, *buffer, dddot, ddddot, u1v2, v1u2;

  char                msg[ 25 ];
  char                msg2[25];
  
  /*
   * External procedures
   * -------------------
   */
  
  extern Integer  mxmynd_(), mxnprc_();
  extern void     mxinit_();

  extern Integer  mapchk_(), count_list();
  extern void     memreq_();
  extern void     pdiff(), xstop_(), pgexit(), mapdif_(), reduce_maps_();

  extern void     pdspevx();
  
  /*
   *  ---------------------------------------------------------------
   *                      Executable Statements
   *  ---------------------------------------------------------------
   */
  
  /*
   *  Get this processor's id, and the number of allocated nodes.
   */
  
  mxinit_();

  me    = mxmynd_();
  nproc = mxnprc_();
  
  strcpy( msg,  "Error in pdspev." );
  
#ifdef DEBUG
   printf("me = %d In pdspev \n", me );
#endif

 /*
  *     Test the input parameters.
  */
  
  linfo = 0;
  
  msize = *n;
  *info = 0;
  
 /*
  *  Quick return if possible.
  */
  
  if ( *n == 0 )
    return;
  
  nvecsZ = count_list( me, mapZ, &msize );
  
  if ( nvecsZ <= 0 )
    return;
  
  for ( k = 0; k < nvecsZ; k++ )
    if ( vecZ[ k ] == NULL )
      *info = -4;
  
   if ( *info != 0 ) {
       linfo = *info;
       fprintf( stderr, " %s me = %d argument %d contains a pointer to NULL. \n",
                msg, me, -linfo);
       xstop_( info );
       return;
   }
  
  /*
   *  ------------------------------------------------
   *  No local errors, compare data across processors.
   *  ------------------------------------------------
   */

  /*
   *  Reduce mapA and mapZ to a single sorted list of processors.
   */

   proclist = iscratch;
   reduce_maps( msize, mapZ, *n, mapZ, 0, mapZ, &nn_proc, proclist );
   
   /* -------------------------------------
    * All input data is good.  Call pdspevx.
    * -------------------------------------
    */
   
   /*
    * Set "extra" arguments so pdspevx  will compute all eigenvalues
    * and eigenvectors.
    */

   m = msize/2;
   /*
     for ( jjj = 1; jjj < m; jjj++ )
     eval[jjj] = eval[2*jjj];
   */
   
   l_good = 1;

   nvecZ = 0;
   cvecZ = 0;
   if ( me == 0 ){
     nvecZ = 1; /* next position */
     cvecZ = 1;
   }
   buffer = scratch;
   for ( k = 1; k < *n; k++ ) { /* reduce *n=2*m vects to m vects */
     if ( mapZ[k] == me ){
       dcopy_(&msize, vecZ[cvecZ], &IONE, buffer, &IONE);
       nvecZ++;
     }
     
     bbcast00((char *) buffer, msize*sizeof(DoublePrecision),
	      k, mapZ[k], nn_proc, proclist);
     
     /*
       compare broadcasted vector with those in good list
       the good vectors overwrites the previous vectors
     */
     
     dddot = 0.;
     
     for ( jjj = 0; jjj < l_good; jjj++ ){
       ccvecZ = 0;
       if ( mapZ[jjj] == me ){
	 /*
	   real part is ortho already; test complex part
	   u = u1 + i u2 = buffer vector
	   v = v1 + i v2 = good vector
	   imag(u.v bar) = -u1v2 + v1u2
	 */
	 
	 u1v2 = ddot_( &m, buffer, &IONE, &vecZ[ccvecZ][m], &IONE);
	 v1u2 = ddot_( &m, &buffer[m], &IONE, &vecZ[ccvecZ][0], &IONE);
	 dddot = max(fabs(-u1v2+v1u2), dddot);
	 ccvecZ++;
       }
     }
     
     gmax00( (char *) &dddot, 1, 5, 16, proclist[0], nn_proc, proclist, &scratch[msize]);
     
     /*
       gmax of dot
     */
     
     if ( ddddot < 1.e-10 ){
       if ( mapZ[l_good] == me ) {
	 eval[cvecZ] = eval[k];
	 dcopy_(&msize, buffer, &IONE, vecZ[cvecZ], &IONE);
	 cvecZ++;
       }
       l_good++;
     }
   }
   
   /*
     good vectors = 
   */
   if ( l_good != m ) {
     printf(" Hell is loose---pdcomplex \n");
     exit(-1);
   }
   
   return;
}



