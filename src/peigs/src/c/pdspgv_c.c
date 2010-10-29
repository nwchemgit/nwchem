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


void pdspgv( ifact, n, vecA, mapA, vecB, mapB, vecZ, mapZ,
	    eval, iscratch, iscsize, dblptr, ibuffsize, scratch, ssize, info)
     Integer            *ifact, *n, *mapA, *mapB, *mapZ, *iscratch, *iscsize, *ibuffsize, *ssize, *info;
     DoublePrecision         **vecA, **vecB, **vecZ, *eval, **dblptr, *scratch;
{
  
  /*
   *
   *  Our parallel version of LAPACK's dspgv.
   *
   *
   *  Purpose
   *  =======
   *
   *  pdspgv computes all of the eigenvalues and eigenvectors of a
   *  real generalized symmetric-definite eigenproblem, of the form
   *
   *  A*x=(lambda)*B*x.
   *
   *  Here A and B are assumed to be symmetric and B is also
   *  positive definite. The matrices A and B use packed storage.
   *
   *  Arguments
   *  =========
   *
   * NOTE: The current code assumes mapA, mapB and mapZ each contain
   *       exactly the same set of processor id's, though not necessarily
   *       in the same order.
   *
   * All arguments are POINTERS to date of the specified type unless
   * otherwise mentioned.  In particular, INTEGER = (Integer *) and
   * DOUBLE PRECISION = (DoublePrecision *)
   *
   *
   * In the following let:
   *    n      = dimension of matrices A and B
   *    me     = this processors id (= mxmynd_())
   *    nprocs = number of allocated processors ( = mxnprc_())
   *    nvecsA = number of entries in mapA equal to me
   *             (= count_list( me, mapA, n ))
   *             number of columns of A this processor has
   *    nvecsB = number of entries in mapB equal to me
   *             (= count_list( me, mapB, n ))
   *             number of vectors this processor has
   *    nvecsZ = number of entries in mapZ equal to me
   *             (= count_list( me, mapZ, n ))
   *
   *-----------------------------------------------------------------
   *
   *  ifact   (input) INTEGER
   *          Specifies whether or not to factor the B matrix.
   *          = 0:  Don't factor B, vecB is already the L in B = L*L'.
   *                (for serial ) or L**(-1) for parallel
   *                which was computed on a previous call to pdspgv.
   *
   *          = 1:  Factor B
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
   *  vecB    (input/output) array of pointers to DoublePrecision (DoublePrecision **)
   *                            dimension ( nvecsB )
   *          On entry, vecB[i], i = 0 to nvecsB-1, points to an
   *          array containing the lower triangular part of the i-th
   *          column of B which is owned by this processor.  The
   *          columns of B owned by this processer are determined by mapB
   *          (See below).
   *
   *          On exit
   *            Let L be the triangular factor L from the Choleski
   *            factorization B = L*L', then
   *            
   *            if (nprocs = 1): vecB contains  L,          same storage as B
   *            else           : vecB contains (L inverse), same storage as B
   *
   *
   *  mapB    (input) INTEGER array, dimension (N)
   *          mapB(i) = the id of the processor which owns column i
   *                    of the B matrix, i = 0 to n-1.
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
   *            The eigenvectors are normalized such that: Z'*B*Z = I.
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
   *  iscratch (workspace) INTEGER array, dimension (??)
   *
   *  iscsize  (input) INTEGER
   *           The number of usable elements in array "iscratch".
   *           Must be >= ???.
   *
   *  dblptr   (workspace) DoublePrecision pointer to DoublePrecision (DoublePrecision **),
   *           dblprt[0] must point to the start of an array consising of
   *           ??? pointers to DoublePrecision.
   *
   *  ibuffsize(input) INTEGER
   *           The number of usable elements in array "dblptr".
   *           Must be >= ???.
   *
   *  scratch  (workspace) DOUBLE PRECISION array, dimension (??)
   *
   *  ssize    (input) INTEGER
   *           The number of usable elements in array "scratch".
   *           Must be >= ???.
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
   *          = 1,             Error choleski factoring B.  In this case
   *                           'choleski' returned a non-zero info, whose
   *                           value was printed to stderr.
   *
   *          = 2,             Error computing inverse of L, the choleski
   *                           factor of B.  In this case
   *                           'inverseL' returned a non-zero info, whose
   *                           value was printed to stderr.
   *
   *          = 3,             Error solving the standard eigenproblem.
   *                           In this case
   *                           'pdspevx' returned a non-zero info, whose
   *                           value was printed to stderr.
   *
   *          Any processor with a negative INFO stops program execution.
   *
   *          If -50 <= INFO < 0, then bad data was passed to this routine
   *                              and no attempt is made to make sure that
   *                              INFO is the same on all processors.
   *        
   *          All other INFO      INFO should be the same on all processors in
   *                              mapA,B,Z.  If this routine calls a routine
   *                              which returns a negative info, then info
   *                              may not be the same on all processors.
   *                              This, however, should never happen.
   *     ---------------------------------------------------------------------
   */
  
  /*
   * Local variables
   * ---------------
   */
  
  Integer             i, j, k, isize, msize, me, nproc;
  Integer             nn_proc, nvecsA, nvecsB, nvecsZ;
  Integer             linfo, maxinfo;
  Integer             irange, ivector, meigval, ilb, iub;
  Integer             *i_scrat, *proclist;
  
  char            msg[ 40 ];
  char            msg2[ 40 ];
  
  DoublePrecision          lb, ub, abstol;
  
/*
 * External procedures
 * -------------------
 */

    extern Integer      mxmynd_(), mxnprc_();

    extern Integer  mapchk_(), count_list();
    extern void     reduce_maps();
    extern void     pdiff(), pgexit();
    extern void     mapdif_();
    extern void     memreq_();
    extern void     pdspgvx();
    extern void     mxinit_(), xstop_();


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
  
  strcpy( msg,  "Error in pdspgv." );
  
#ifdef DEBUG1
   fprintf(stderr, "me = %d In pdspgv \n", me );
#endif

  /*
   *     Test the input parameters.
   */
  
    linfo = 0;
    
    if ( ifact  == NULL )
      linfo = -1;
    else if ( n == NULL )
      linfo = -2;
    else if ( vecA == NULL )
      linfo = -3;
    else if ( mapA == NULL )
      linfo = -4;
   else if ( vecB == NULL )
      linfo = -5;
   else if ( mapB == NULL )
      linfo = -6;
   else if ( vecZ == NULL )
      linfo = -7;
   else if ( mapZ == NULL )
      linfo = -8;
   else if ( eval == NULL )
      linfo = -9;
   else if ( iscratch == NULL )
     linfo = -10;
   else if ( iscsize == NULL )
      linfo = -11;
   else if ( dblptr == NULL )
      linfo = -12;
   else if ( ibuffsize == NULL )
      linfo = -13;
   else if ( scratch == NULL )
      linfo = -14;
   else if ( ssize == NULL )
     linfo = -15;
   else if ( info == NULL )
     linfo = -16;
    
   if ( linfo != 0 ){
      if ( info != NULL )
	*info = linfo;
      
      fprintf( stderr, " %s me = %d argument %d is a pointer to NULL. \n",
               msg, me, -linfo );
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

    if ( *ifact < 0 || *ifact > 1 )
      *info = -1;
    else if ( *n < 0 )
      *info = -2;
    else if( mapchk_( n, mapA ) != 0 )
      *info = -4;
    else if( mapchk_( n, mapB ) != 0 )
      *info = -6;
    else if( mapchk_( n, mapZ ) != 0 )
      *info = -8;
    
    if ( *info != 0 ) {
        linfo = *info;
        fprintf( stderr, " %s me = %d argument %d has an illegal value. \n",
                 msg, me, -linfo);
        xstop_( info );
	return;
    }
    
    if ( *info == 0 ) {
	/*
	 * Count the number of columns of A, B and Z owned by this processor.
	 * Must own something.
	 */
	
	nvecsA = count_list( me, mapA, &msize );
	nvecsB = count_list( me, mapB, &msize );
	nvecsZ = count_list( me, mapZ, &msize );
	
	if ( nvecsA + nvecsB + nvecsZ <= 0 )
	  return;
	
        for ( k = 0; k < nvecsZ; k++ )
	  if ( vecZ[ k ] == NULL )
	    *info = -7;
	
	for ( k = 0; k < nvecsB; k++ )
	  if ( vecB[ k ] == NULL )
	    *info = -5;
	
	for ( k = 0; k < nvecsA; k++ )
	  if ( vecA[ k ] == NULL )
	    *info = -3;
    }
    
    if ( *info != 0 ) {
        linfo = *info;
        fprintf( stderr, " %s me = %d argument %d contains a pointer to NULL. \n",
                 msg, me, -linfo);
        xstop_( info );
	return;
    }
    
    /*
     *  For now REQUIRE mapA, mapB and mapZ to contain exactly the same
     *  set of processors, but not necessarily in the same order.
     */
    
    mapdif_( n, mapA, mapB, iscratch, &linfo );
    if ( linfo != 0 )
      *info = -6;

    mapdif_( n, mapA, mapZ, iscratch, &linfo );
    if ( linfo != 0 )
      *info = -8;
    
    if ( *info != 0 ) {
        fprintf( stderr, " %s me = %d mapA,B,Z differ \n", msg, me );
        xstop_( info );
	return;
    }
    
    /*
     * check memory allocation
     * i is the integer scratch size
     * j is the DoublePrecision precision scratch size
     * k is the pointer scratch size
     */

    irange = 0;
  
    memreq_( &irange, n, mapA, mapB, mapZ, &i, &j, &k, iscratch );
  
    linfo = 0;
    if ( *iscsize < i )
      linfo = -11;
    else if ( *ibuffsize < k )
      linfo = -13;
    else if ( *ssize < j )
      linfo = -15;
  
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
     *  Reduce mapA, mapB and mapZ to a single sorted list of processors.
     */

#ifdef DEBUG7
   printf("me = %d pdspgv 1 \n", me );
#endif
    
    proclist = iscratch;
    reduce_maps( *n, mapA, *n, mapB, *n, mapZ, &nn_proc, proclist );

    i_scrat = iscratch + nn_proc;

    /*
     *  Check scaler inputs.
     */
    
    i_scrat[ 0 ] = *n;
    i_scrat[ 1 ] = *ifact;
    
  isize = 2 * sizeof( Integer );
  strcpy( msg2,  "n or ifact " );
  pdiff( &isize, (char *) i_scrat, proclist, &nn_proc, i_scrat+2, msg, msg2, &linfo );

  pgexit( &linfo, msg, proclist, &nn_proc, scratch );
  
  if ( linfo != 0 ){
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

  strcpy( msg2,  "mapB " );
  pdiff( &isize, (char *) mapB, proclist, &nn_proc, i_scrat, msg, msg2, &linfo );
  maxinfo = max( maxinfo, linfo );

  strcpy( msg2,  "mapZ " );  
  pdiff( &isize, (char *) mapZ, proclist, &nn_proc, i_scrat, msg, msg2, &linfo );
  maxinfo = max( maxinfo, linfo );

  linfo = maxinfo;

  pgexit( &linfo, msg, proclist, &nn_proc, scratch );

  if ( linfo != 0 ){
    *info = -51;
    return;
  }
  
  /* --------------------------------------
   * All input data is good.  Call pdspgvx.
   * --------------------------------------
   */
  
  /*
   * Set "extra" arguments so pdspgvx  will compute all eigenvalues
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

#ifdef DEBUG7
   printf("me = %d pdspgv 2\n", me );
#endif
  
  pdspgvx( ifact, &ivector, &irange, n, vecA, mapA, vecB, mapB,
          &lb, &ub, &ilb, &iub, &abstol, &meigval,
          vecZ, mapZ, eval, iscratch, iscsize,
          dblptr, ibuffsize, scratch, ssize, info);
  
#ifdef DEBUG7
   printf("me = %d Exiting pdspgv \n", me );
#endif

  return;
}

