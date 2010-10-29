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
#include <stdlib.h>
#include <memory.h>
#include <time.h>

#include "globalp.c.h"

#ifdef TIMING
#include "timing.h"
#endif

#define max(a,b) ((a) > (b) ? (a) : (b))
#define min(a,b) ((a) < (b) ? (a) : (b))

void pdspgvx( ifact, ivector, irange, n, vecA, mapA, vecB, mapB,
	     lb, ub, ilb, iub, abstol, meigval, vecZ, mapZ, eval,
	     iscratch, iscsize, dblptr, ibuffsize, scratch, ssize, info)
     
     Integer            *ifact, *ivector, *irange, *n, *mapA, *mapB, *ilb, *iub,
       *meigval, *mapZ, *iscratch, *iscsize,
       *ibuffsize, *ssize, *info;
     
     DoublePrecision         *lb, *ub, *abstol, *eval, *scratch;
     DoublePrecision         **vecA, **vecB, **vecZ, **dblptr;

{

/*
 *
 *  A parallel version of LAPACK's dspgv, but extended to allow for
 *  the computation of selected eigenvalues and, optionally, eigenvectors.
 *
 *
 *  Purpose
 *  =======
 *
 *  pdspgvx computes some or all of the eigenvalues and, optionally,
 *  eigenvectors of a real generalized symmetric-definite eigenproblem,
 *  of the form
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
 *                which was computed on a previous call to pdspgvx.
 *
 *          = 1:  Factor B
 *
 *  ivector (input) INTEGER
 *          = 0:    Compute eigenvalues only;
 *          = 1:    Compute eigenvalues and eigenvectors.
 *
 *  irange  (input) INTEGER
 *          = 1:   all eigenvalues will be found;
 *          = 2:   all eigenvalues in the half-open interval (lb, ub]
 *                 will be found;
 *          = 3:   the ilb-th through iub-th eigenvalues will be found.
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
 *            Let L be the triangular factor L from the Cholesky
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
 *  lb      (input) DOUBLE PRECISION
 *          If IRANGE=2,  the lower bound of the interval to be searched
 *          for eigenvalues.  Not referenced if IRANGE = 1 or 3, but
 *          must not be a pointer to NULL.
 *
 *  ub      (input) DOUBLE PRECISION
 *          If IRANGE=2,  the upper bound of the interval to be searched
 *          for eigenvalues.  Not referenced if IRANGE = 1 or 3, but
 *          must not be a pointer to NULL.
 *
 *  ilb     (input) INTEGER
 *          If IRANGE=3,  the index (from smallest to largest) of the
 *          smallest eigenvalue to be returned.  ilb >= 1.
 *          Not referenced if IRANGE = 1 or 2, but
 *          must not be a pointer to NULL.
 *
 *  iub     (input) INTEGER
 *          If RANGE=3, the index (from smallest to largest) of the
 *          largest eigenvalue to be returned.  min(ilb,N) <= iub <= N.
 *          Not referenced if IRANGE = 1 or 2, but
 *          must not be a pointer to NULL.
 *
 *  abstol  (input) DOUBLE PRECISION
 *          The absolute error tolerance for the eigenvalues.
 *          An approximate eigenvalue is accepted as converged
 *          when it is determined to lie in an interval [a,b]
 *          of width less than or equal to
 *
 *                  ABSTOL + EPS *   max( |a|,|b| ) ,
 *
 *          where EPS is the machine precision.  If ABSTOL is less than
 *          or equal to zero, then  EPS*|T|  will be used in its place,
 *          where |T| is the 1-norm of the tridiagonal matrix obtained
 *          by reducing matrix A to tridiagonal form.
 *
 *          See "Computing Small Singular Values of Bidiagonal Matrices
 *          with Guaranteed High Relative Accuracy," by Demmel and
 *          Kahan, LAPACK Working Note #3.
 *
 *  meigval (output) INTEGER
 *          The total number of eigenvalues found.  0 <= meigval <= N.
 *          If IRANGE = 1,   M = N, and if RANGE = 3,   M = IUB-ILB+1.
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
 *---------------------------------------------------------------------
 */

/*
 * Local variables
 * ---------------
 */

  static Integer G_TYPE = 1111;
  
  Integer             k, isize, msize, mapZ0, me, nproc,
                      nn_proc, neigval, nvecsA, nvecsB, nvecsZ,
                      linfo, maxinfo, kinverse, i, j, itype;
  
  Integer            *proclist, *mapvA, *mapvB, *i_scrat;

  Integer            *g_proclist, n_glob_proc, num_procs;

  
  
  char            msg[ 40 ];
  char            msg2[ 40 ];
  
  DoublePrecision         *d_scrat, *dtmp_ptr;
  DoublePrecision        **buff_ptr;
  
  /*
   * External procedures
   * -------------------
   */
  
  extern Integer  mxmynd_(), mxnprc_();
  extern Integer  mapchk_(), count_list(), fil_mapvec_(), reduce_list2();

  extern void     reduce_maps();
  extern void     xstop_(), pdiff();
  extern void     pgexit();
  extern void     mapdif_();
  extern void     choleski(), mxm5x(), inverseL();
  extern void     forwardLU_(), forwardLL_(), pdspevx(), upperUF_(),
                  lsl_conjugation2(), mxinit_();
  extern void     mdiff2_(), mdif2b_();
  extern void     memreq_();
  extern void     bbcast00();
  extern void     gi_sum();
  
#ifdef TIMING
  extern TIMINGG test_timing;
  DoublePrecision mxclock_();
  DoublePrecision  t1, t2, djunk;
#endif
  
  
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
  strcpy( msg,  "Error in pdspgvx." );

#ifdef DEBUG1
      fprintf(stderr, "me = %d In pdspgvx. \n", me );
#endif
#ifdef TIMING
  mxclock_();
#endif
  
  /*
   *     Test the input parameters.
   */
  
  linfo = 0;
  
  if ( ifact  == NULL )
    linfo = -1;
  else if ( ivector == NULL )
    linfo = -2;
  else if ( irange == NULL )
    linfo = -3;
  else if ( n == NULL )
      linfo = -4;
  else if ( vecA == NULL )
    linfo = -5;
  else if ( mapA == NULL )
    linfo = -6;
  else if ( vecB == NULL )
    linfo = -7;
  else if ( mapB == NULL )
    linfo = -8;
  else if ( lb == NULL )
    linfo = -9;
  else if ( ub == NULL )
    linfo = -10;
  else if ( ilb == NULL )
    linfo = -11;
  else if ( iub == NULL )
    linfo = -12;
  else if ( abstol == NULL )
    linfo = -13;
  else if ( meigval == NULL )
    linfo = -14;
  else if ( vecZ == NULL )
    linfo = -15;
  else if ( mapZ == NULL )
    linfo = -16;
  else if ( eval == NULL )
    linfo = -17;
  else if ( iscratch == NULL )
    linfo = -18;
  else if ( iscsize == NULL )
    linfo = -19;
  else if ( dblptr == NULL )
    linfo = -20;
  else if ( ibuffsize == NULL )
    linfo = -21;
  else if ( scratch == NULL )
    linfo = -22;
  else if ( ssize == NULL )
    linfo = -23;
  else if ( info == NULL )
    linfo = -24;
  
  if ( linfo != 0 ){
    if ( info != NULL )
      *info = linfo;
      
    fprintf( stderr, " %s me = %d argument %d is a pointer to NULL. \n",
             msg, me, -linfo );
    xstop_( &linfo );
    return;
  }
  
  msize    = *n;
  *info    = 0;
  *meigval = 0;
  
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
  else if ( *ivector < 0  || *ivector > 1 )
    *info = -2;
  else if ( *irange < 1  || *irange > 3 )
    *info = -3;
  else if ( *n < 0 )
    *info = -4;
  else if ( *irange == 2  && *lb >= *ub )
    *info = -10;
  else if ( *irange == 3  && *ilb < 1 )
    *info = -11;
  else if ( *irange == 3  && ( *iub < *ilb  ||  *iub > *n ) )
    *info = -12;
  else if( mapchk_( n, mapA ) != 0 )
    *info = -6;
  else if( mapchk_( n, mapB ) != 0 )
    *info = -8;
  else if( mapchk_( n, mapZ ) != 0 )
    *info = -16;
    
  if ( *info != 0 ) {
      linfo = *info;
      fprintf( stderr, " %s me = %d argument %d has an illegal value. \n",
               msg, me, -linfo);
      xstop_( info );
      return;
  }
    
   /*
    * Count the number of columns of A, B and Z owned by this processor.
    */
    
    nvecsA = count_list( me, mapA, &msize );
    nvecsB = count_list( me, mapB, &msize );
    nvecsZ = count_list( me, mapZ, &msize );
    
    if ( nvecsA + nvecsB + nvecsZ <= 0 )
      return;
    
    if ( *info == 0 ) {
      dtmp_ptr = vecZ[0];
      for ( k = 0; k < nvecsZ; k++ )
	if ( dtmp_ptr++ == NULL ) {
	  *info = -15;
	  break;
	}
    }
    
    if ( *info == 0 ) {
      dtmp_ptr = vecB[0];
      for ( k = 0; k < nvecsB; k++ )
	if ( dtmp_ptr++ == NULL ) {
	  *info = -7;
	  break;
	}
    }
      
    if ( *info == 0 ) {
      dtmp_ptr = vecA[0];
      for ( k = 0; k < nvecsA; k++ )
	if ( dtmp_ptr++ == NULL ) {
	  *info = -5;
	  break;
	}
    }
  
  if ( *info != 0 ) {
      linfo = *info;
      fprintf( stderr, " %s me = %d argument %d contains a pointer to NULL. \n",
               msg, me, -linfo);
      xstop_( info );
      return;
  }
    
  /*
   *  REQUIRE mapA, mapB and mapZ to contain exactly the same
   *  set of processors, but not necessarily in the same order.
   */
    
  mapdif_( n, mapA, mapZ, iscratch, &linfo );
  if ( linfo != 0 )
    *info = -16;
    
  mapdif_( n, mapA, mapB, iscratch, &linfo );
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
  
  itype = 0;
  memreq_( &itype, n, mapA, mapB, mapZ, &i, &j, &k, iscratch );
  
  linfo = 0;
  if ( *iscsize < i )
    linfo = -19;
  else if ( *ibuffsize < k )
    linfo = -21;
  else if ( *ssize < j )
    linfo = -23;
  
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
  
  proclist = iscratch;

  reduce_maps( *n, mapA, *n, mapB, *n, mapZ, &nn_proc, proclist );

  i_scrat = iscratch + nn_proc;

  /*
   *  Check Integer scaler inputs.
   */
  
  i_scrat[ 0 ] = *ifact;
  i_scrat[ 1 ] = *ivector;
  i_scrat[ 2 ] = *irange;
  i_scrat[ 3 ] = *n;
  i_scrat[ 4 ] = *ilb;
  i_scrat[ 5 ] = *iub;
  
  isize = 6 * sizeof( Integer );
  strcpy( msg2,  "ifact,ivector,irange,n,ilb,or iub " );
  pdiff( &isize, (char *) i_scrat, proclist, &nn_proc, i_scrat+6, msg, msg2, &linfo );
  
  pgexit( &linfo, msg, proclist, &nn_proc, scratch );
  
  if ( linfo != 0 ){
    *info = -51;
    return;
  }
  
  /*
   *  Check remaining inputs.
   */
  
  maxinfo = 0;
  
  isize = msize * sizeof( Integer );
  
  strcpy( msg2,  "mapA " );
  pdiff( &isize, (char *) mapA, proclist, &nn_proc, i_scrat, msg, msg2, &linfo );
  maxinfo = max( maxinfo, linfo );
  
  strcpy( msg2,  "mapB " );
  pdiff( &isize, (char *) mapB, proclist, &nn_proc, i_scrat, msg, msg2, &linfo );
  maxinfo = max( maxinfo, linfo );
  
  strcpy( msg2,  "mapZ " );
  pdiff( &isize, (char *) mapZ, proclist, &nn_proc, i_scrat, msg, msg2, &linfo );
  maxinfo = max( maxinfo, linfo );
  
  isize        = 3 * sizeof( DoublePrecision );
  scratch[ 0 ] = *lb;
  scratch[ 1 ] = *ub;
  scratch[ 2 ] = *abstol;
  
  strcpy( msg2,  "lb,ub,or abstol " );
  pdiff( &isize, (char *) scratch, proclist, &nn_proc, (Integer *) (scratch+3), msg, msg2, &linfo );
  maxinfo = max( maxinfo, linfo );
  
  linfo = maxinfo;
  
  pgexit( &linfo, msg, proclist, &nn_proc, scratch );
  
  if ( linfo != 0 ){
    *info = -51;
    return;
  }
  
  /* ----------------------------------------
   * All input data is good. Start computing. 
   * ----------------------------------------
   */
  
  /*
   *  Use L inverse when number of allocated processors is greater than 1.
   */
  
  /*
   *  Any modifications to the decision to use L inverse must be made in
   *  memreq.  memreq.c currently always uses L inverse
   */
  
  if ( nproc == 1 )
    kinverse = 0;
  else
    kinverse = 1;

  kinverse = 1;
  
#ifdef TIMING
  kinverse = 1;
#endif
  
    /*
     * Initialize workspace.
     */

    /* i_scrat set above */

    d_scrat  = scratch;
    buff_ptr = dblptr;

    n_glob_proc = nn_proc;
    g_proclist = proclist;

    proclist = NULL;
  
    mapvB = i_scrat;
    i_scrat += nvecsB;

    mapvA = i_scrat;
    i_scrat += nvecsA;

    proclist = i_scrat;
    i_scrat += nproc;
    
    /*
     * Fill in mapv* arrays.
     */
  
  fil_mapvec_(&me, &msize, mapA, mapvA);
  fil_mapvec_(&me, &msize, mapB, mapvB);

  mapZ0 = *mapZ;    
  
  if (*ifact == 1){
    if ( nvecsB > 0 ){
      
      /*
       *  Factor B
       *  --------
       */
      
#ifdef TIMING
      t1 = mxclock_();
#endif
      
      choleski( &msize, vecB, mapB, i_scrat, d_scrat, info );
      
      if( *info != 0 ) {
        fprintf(stderr, " %s me = %d choleski returned info = %d \n", msg, me, *info );
        *info = 1;
      }
      
#ifdef TIMING
      mxsync_();
      t2 = mxclock_();
      test_timing.choleski = t2 - t1;
#endif
      
      /*
       *  Invert L
       *  --------
       */
      
      if ( kinverse == 1  &&  *info == 0 ){
	
#ifdef TIMING
	t1 = mxclock_();
#endif
	
	inverseL( &msize,  vecB, mapB, i_scrat, d_scrat, info);
	
#ifdef TIMING
	mxsync_();
	t2 = mxclock_();
	test_timing.inverse = t2 - t1;
#endif

        if( *info != 0 ) {
          fprintf(stderr, " %s me = %d inverseL returned info = %d \n", msg, me, *info );
          *info = 2;
        }
      
      }
    }
    
    /*
     * Send info to all processors in {mapA + mapZ} - {map B}
     */
    
    mdiff2_( n, mapA, n, mapZ, n, mapB, i_scrat, &num_procs ); 

    if( num_procs > 0 ){
      i_scrat[num_procs] = mapB[0];
      num_procs++;
      bbcast00((char *) info, sizeof(Integer), 2, mapB[0], num_procs, i_scrat);
    }
    
    if ( *info != 0 )
      return;
    
  }
  
  if ( nvecsA + nvecsB > 0  &&  kinverse == 0 ){
    forwardLU_( &msize, mapB, mapvB, vecB, mapA, mapvA, vecA, d_scrat,
	       &nn_proc, proclist, i_scrat, d_scrat + 2 * msize,
	       buff_ptr);
    
    forwardLL_( &msize, mapB, mapvB, vecB, mapA, mapvA, vecA,
	       d_scrat, &nn_proc, proclist, i_scrat);
    
  }
  else if ( nvecsA + nvecsB > 0  &&  kinverse == 1 ){

#ifdef TIMING
    t1 = mxclock_();
#endif
    
    lsl_conjugation2( &msize, vecA, mapA, vecB, mapB, i_scrat, d_scrat,buff_ptr );
    
#ifdef TIMING
    mxsync_();
    t2 = mxclock_();
    test_timing.conjug = t2 - t1;
#endif

  }
  
  /*
   * Set num_procs  = (# of processors in mapB, but not mapA
   *                  or mapZ) + 1.
   *
   *     proclist[0:num_procs-1] = list of procssors in
   *                               (mapB - {mapA + mapZ}) U {mapZ[0]}.
   *
   * This data is for bbcast0 below.  We must do this before calling
   * pdspevx since that routine modifies part of mapZ.
   */
  
   mdif2b_( n, mapB, n, mapA, n, mapZ, proclist, &num_procs ); 

   if( num_procs > 0 ){
     proclist[num_procs] = mapZ0;
     num_procs++;
   }
  
  /*
   * iscsize, ibuffsize, ssize need to be modified.
   */
  
  if (nvecsA + nvecsZ > 0) {

#ifdef TIMING
    t1 = mxclock_();
#endif
    
    pdspevx( ivector, irange, n, vecA, mapA, lb, ub, ilb, iub, abstol,
	    meigval, vecZ, mapZ, eval, i_scrat, iscsize,
	    buff_ptr, ibuffsize, d_scrat, ssize, info);
    
#ifdef TIMING
    mxsync_();
    t2 = mxclock_();
    test_timing.pdspevx = t2 - t1;
#endif

    if( *info != 0 ) {
      fprintf(stderr, " %s me = %d pdspevx returned info = %d \n", msg, me, *info );
      *info = 3;
	exit(-1);
    }
  }

  /*
     out global list pdspgvx sync
     */  
  
  i = 0;
  gi_sum( &i, 1, G_TYPE,  0, n_glob_proc, g_proclist, d_scrat);
  
  /*
   * Send info, mapZ, meigval and eval to all processors in
   * {mapB} - {mapA + mapZ}
   */
  
  if (num_procs > 0){
    bbcast00((char *) info, sizeof(Integer), 2, mapZ0, num_procs, proclist);
    bbcast00((char *) mapZ, msize*sizeof(Integer), 3, mapZ0, num_procs, proclist);
    bbcast00((char *) meigval, sizeof(Integer), 4, mapZ0, num_procs, proclist);
    bbcast00((char *) eval, *meigval* sizeof(DoublePrecision), 5, mapZ0, num_procs, proclist);
  }
  
  i = 0;
  gi_sum( &i, 1, G_TYPE,  0, n_glob_proc, g_proclist, d_scrat);
  
  neigval = *meigval;
  
  nvecsZ = count_list( me, mapZ, &msize );
  
  if ( *info != 0 )
    return;
  
  /*
   * Exit if no eigenvectors need to be back-transformed, otherwise,
   * back-transform the eigenvectors.
   */
  
  if ( neigval == 0  ||  *ivector == 0 )
    return;
  
  if (nvecsB + nvecsZ > 0) {
      
#ifdef TIMING
      t1 = mxclock_();
#endif
      if ( kinverse == 0 )
	upperUF_( &msize, vecB, mapB, &neigval, vecZ, mapZ, i_scrat, d_scrat);
      else
	mxm5x( &msize, vecB, mapB, &neigval, vecZ, mapZ, i_scrat, d_scrat );
      
#ifdef TIMING
      mxsync_();
      t2 = mxclock_();
      test_timing.mxm5x = t2 - t1;
#endif

  }

#ifdef DEBUG1
      fprintf(stderr, "me = %d Exiting pdspgvx. \n", me );
#endif

/*
	for ( i = 0; i < *n; i++ )
	       printf(" eval i %d eval %g \n", i, eval[i]);
*/


  return;
}
