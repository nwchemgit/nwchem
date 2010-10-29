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

void pdspev(n, vecA, mapA, vecZ, mapZ, eval, iscratch, iscsize,
            dblptr, ibuffsize, scratch, ssize, info)
     
     Integer            *n, *mapA, *mapZ, *iscratch, *iscsize, *ibuffsize, *ssize, *info;
     DoublePrecision         **vecA, **vecZ, *eval, **dblptr, *scratch;
{
  
/*
 *
 *  Our parallel version of LAPACK's dspev.
 *
 *
 *  Purpose
 *  =======
 *
 *  pdspev computes all of the eigenvalues and eigenvectors of a
 *  real symmetric eigenproblem, of the form
 *
 *  A*x=(lambda) * x.
 *
 *  Here A is assumed to be symmetric.
 *  A uses packed storage.
 *
 *  Arguments
 *  =========
 *
 * The current code assumes mapA and mapZ each contain
 * exactly the same set of processor id's, though not necessarily
 * in the same order.
 *
 * All arguments are POINTERS to date of the specified type unless
 * otherwise mentioned.  In particular, INTEGER = (Integer *) and
 * DOUBLE PRECISION = (DoublePrecision *)
 *
 *
 * In the following let:
 *    n      = dimension of matrix A
 *    me     = this processors id (= mxmynd_())
 *    nprocs = number of allocated processors ( = mxnprc_())
 *    nvecsA = number of entries in mapA equal to me
 *             (= count_list( me, mapA, n ))
 *    nvecsZ = number of entries in mapZ equal to me
 *             (= count_list( me, mapZ, n ))
 *
 *-----------------------------------------------------------------
 *
 *  n       (input) INTEGER
 *          The number of rows and columns of the matrices A and B.
 *          N >= 0.
 *
 *  vecA    (input/workspace) array of pointers to DoublePrecision (DoublePrecision **)
 *                            dimension ( nvecsA )
 *          On entry, vecA[i], i = 0 to nvecsA-1, points to an
 *          array containing the lower triangular part of the i-th
 *          column of A which is owned by this processor.  The
 *          columns of A owned by this processer are determined by mapA
 *          (See below).
 *
 *          On exit, the contents of matrixA are destroyed.
 *
 *  mapA    (input) INTEGER array, dimension (N)
 *          mapA(i) = the id of the processor which owns column i
 *                    of the A matrix, i = 0 to n-1.
 *
 *  vecZ    (output) array of pointers to DoublePrecision (DoublePrecision **)
 *                   dimension ( nvecsZ )
 *          On entry, vecZ[i], i = 0 to nvecsZ-1, should point to an
 *          array of length n.
 *
 *          On exit:
 *
 *            vecZ[i], i = 0 to nvecsZ-1, points to the i-th eigenvector
 *            (as determined by the exit values in mapZ) owned by this
 *            processor.
 *
 *            The eigenvectors are normalized such that: Z'*Z = I.
 *
 *  mapZ    (input/output) INTEGER array, dimension (N)
 *          On entry:
 *
 *          mapZ(i) = the id of a processor which has room for the
 *                    i-th eigenvector, i = 0 to n-1.
 *          On exit:
 *            mapZ(i) = the id of the processor which actually owns the i-th
 *                      eigenvector, i = 0 to n-1.
 *
 *  eval    (output) DOUBLE PRECISION array, dimension (N)
 *          If INFO = 0, the eigenvalues  of the matrix
 *          in no particular order.
 *
 *  iscratch (workspace) INTEGER array, dimension
 *           Must be >= that returned from utility routine memreq_
 *
 *  iscsize  (input) INTEGER
 *           The number of usable elements in array "iscratch".
 *           Must be >= that returned from utility routine memreq_
 *
 *  dblptr   (workspace) DoublePrecision pointer to DoublePrecision (DoublePrecision **),
 *           dblprt[0] must point to the start of an array consising of
 *           pointers to DoublePrecision (DoublePrecision *).
 *           Must be >= that returned from utility routine memreq_
 *
 *  ibuffsize(input) INTEGER
 *           The number of usable elements in array "dblptr".
 *           Must be >= that returned from utility routine memreq_
 *
 *  scratch  (workspace) DOUBLE PRECISION array, dimension
 *           Must be >= that returned from utility routine memreq_
 *
 *  ssize    (input) INTEGER
 *           The number of usable elements in array "scratch".
 *           Must be >= that returned in utility routine memreq_
 *
 *  INFO    (output) INTEGER
 *
 *          = 0:  successful exit.
 *
 *          -50 <= INFO < 0, the -INFOth argument had an illegal value.
 *
 *          INFO = -51,      Input data which must be the same on all
 *                           processors is NOT the same on all processors
 *                           in mapZ.
 *
 *          = 1,             Error computing the eigenvalues of the
 *                           tridiagonal eigenproblem.  In this case
 *                           'pstebz_' returned a non-zero info, whose
 *                           value was printed to stderr.
 *
 *          = 2,             Error computing the eigenvectors of the
 *                           tridiagonal eigenproblem.  In this case
 *                           'pstein' returned a non-zero info, whose
 *                           value was printed to stderr.
 *
 *          Any processor with a negative INFO stops program execution.
 *
 *          If -50 <= INFO < 0, then bad data was passed to this routine
 *                              and no attempt is made to make sure that
 *                              INFO is the same on all processors.
 *        
 *          All other INFO      INFO should be the same on all processors in
 *                              mapA,Z.  If this routine calls a routine
 *                              which returns a negative info, then info
 *                              may not be the same on all processors.
 *                              This, however, should never happen.
 *---------------------------------------------------------------------
 */
  
  /*
   * Local variables
   * ---------------
   */
  Integer             k, me, nn_proc, msize, irange, ivector, meigval, ilb,
                      iub, nproc, isize, nvecsA, nvecsZ, linfo, maxinfo,
                      i, j, itmp;
  Integer             *i_scrat, *proclist;

  DoublePrecision     lb, ub, abstol;

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
  extern void     pdiff(), xstop_(), pgexit(), mapdif_(), reduce_maps();

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
  
  if ( n == NULL )
    linfo = -1;
  else if ( vecA == NULL )
    linfo = -2;
  else if ( mapA == NULL )
    linfo = -3;
  else if ( vecZ == NULL )
    linfo = -4;
  else if ( mapZ == NULL )
    linfo = -5;
  else if ( eval == NULL )
    linfo = -6;
  else if ( iscratch == NULL )
    linfo = -7;
  else if ( iscsize == NULL )
    linfo = -8;
  else if ( dblptr == NULL )
    linfo = -9;
  else if ( ibuffsize == NULL )
    linfo = -10;
  else if ( scratch == NULL )
    linfo = -11;
  else if ( ssize == NULL )
    linfo = -12;
  else if ( info == NULL )
    linfo = -13;
  
  if ( linfo != 0 ){
    if ( info != NULL )
      *info = linfo;
    
    fprintf( stderr, " %s me = %d argument %d is a pointer to NULL. \n", msg, me, -linfo );
    xstop_( &linfo );
    return;
  }
  
  msize = *n;
  *info = 0;
  
 /*
  *  Quick return if possible.
  */
  
  if ( *n == 0 )
    return;
  
 /*
  *  Continue error checking.
  */
  
  if ( *n < 0 )
    *info = -1;
  else if( mapchk_( n, mapA ) != 0 )
    *info = -3;
  else if( mapchk_( n, mapZ ) != 0 )
    *info = -5;
    
  if ( *info != 0 ) {
      linfo = *info;
      fprintf( stderr, " %s me = %d argument %d has an illegal value. \n",
               msg, me, -linfo);
      xstop_( info );
      return;
  }

  /*
   * Count the number of columns of A and Z owned by this processor.
   */
    
   nvecsA = count_list( me, mapA, &msize );
   nvecsZ = count_list( me, mapZ, &msize );
    
   if ( nvecsA + nvecsZ <= 0 )
     return;

   for ( k = 0; k < nvecsZ; k++ )
     if ( vecZ[ k ] == NULL )
       *info = -4;
	
   for ( k = 0; k < nvecsA; k++ )
     if ( vecA[ k ] == NULL )
       *info = -2;

   if ( *info != 0 ) {
       linfo = *info;
       fprintf( stderr, " %s me = %d argument %d contains a pointer to NULL. \n",
                msg, me, -linfo);
       xstop_( info );
       return;
   }
  
  /*
   *  REQUIRE mapA, mapZ to contain exactly the same
   *  set of processors, but not necessarily in the same order.
   */
    
   mapdif_( n, mapA, mapZ, iscratch, &linfo );
   if ( linfo != 0 )
     *info = -5;
    
   if ( *info != 0 ) {
       fprintf( stderr, " %s me = %d mapA,Z differ \n", msg, me );
       xstop_( info );
       return;
   }
    
   /*
    * check memory allocation
    * i is the integer scratch size
    * j is the DoublePrecision precision scratch size
    * k is the pointer scratch size
    */
    
   itmp = 1;
   memreq_( &itmp, n, mapA, mapZ, mapZ, &i, &j, &k, iscratch );
  
   linfo = 0;
   if ( *iscsize < i )
     linfo = -8;
   else if ( *ibuffsize < k )
     linfo = -10;
   else if ( *ssize < j )
     linfo = -12;
   
   if ( linfo != 0 ) {
     *info = linfo;
     fprintf( stderr, " %s me = %d Insufficient workspace provided. \n", msg, me );
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
   reduce_maps( *n, mapA, *n, mapZ, 0, mapZ, &nn_proc, proclist );

   i_scrat = iscratch + nn_proc;

  /*
   *  Check scaler inputs.
   */

   i_scrat[ 0 ] = *n;

   isize = 1 * sizeof( Integer );
   strcpy( msg2,  "n " );
   pdiff( &isize, (char *) i_scrat, proclist, &nn_proc, i_scrat+2, msg, msg2 , &linfo );
   
   pgexit( &linfo, msg, proclist, &nn_proc, scratch );
   
   if ( linfo != 0 ) {
     *info = -51;
     return;
   }
   
   /*
    *  Check rest of inputs.
    */
   
   maxinfo = 0;
   
   isize   = msize * sizeof( Integer );
   strcpy( msg2,  "mapA " );
   pdiff( &isize, (char *) mapA, proclist, &nn_proc, i_scrat, msg, msg2, &linfo );
   maxinfo = max( maxinfo, linfo );
   
   strcpy( msg2,  "mapZ " );
   pdiff( &isize, (char *) mapZ, proclist, &nn_proc, i_scrat, msg, msg2, &linfo );
   maxinfo = max( maxinfo, linfo );
   
   linfo = maxinfo;
   
   pgexit( &linfo, msg, proclist, &nn_proc, scratch );
   
   if ( linfo != 0 ) {
     *info = -51;
     return;
   }
   
   
   /* -------------------------------------
    * All input data is good.  Call pdspevx.
    * -------------------------------------
    */
   
   /*
    * Set "extra" arguments so pdspevx  will compute all eigenvalues
    * and eigenvectors.
    */
  
  irange  =  1;
  ivector =  1;
  meigval =  0;
  ilb     =  0;
  iub     = -1;
  lb      = (DoublePrecision )  0.0e0;
  ub      = (DoublePrecision ) -1.0e0;
  abstol  = (DoublePrecision )  0.0e0;
  
  
  /*
    lll = 0;
    for ( k = 0; k < *n; k++ ) {
    if ( mapZ[k] == me ) {
    for ( jjj = 0; jjj < *n; jjj++ )
    vecZ[lll][jjj] = 0.0e0;
    lll++;
    }
    }
  */
  
  pdspevx( &ivector, &irange, n, vecA, mapA, &lb, &ub, &ilb, &iub, &abstol,
	   &meigval, vecZ, mapZ, eval, iscratch, iscsize,
	   dblptr, ibuffsize, scratch, ssize, info);
  
#ifdef DEBUG99
  printf("me = %d Exiting pdspev \n", me );

  lll = 0;
  for ( k = 0; k < *n; k++ ) {
    if ( mapZ[k] == me ) {
      for ( jjj = 0; jjj < *n; jjj++ )
	printf(" %d  %d %g \n", lll, jjj, vecZ[lll][jjj]);
      lll++;
    }
  }

  for ( jjj = 0; jjj < *n; jjj++ )
    printf(" eval %d %g \n", jjj, eval[jjj]);

#endif
  
  

    return;
}



