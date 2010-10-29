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
#include <time.h>

#include "globalp.c.h"

#ifdef TIMING
#include "timing.h"
#endif

#define max(a,b) ((a) > (b) ? (a) : (b))
#define min(a,b) ((a) < (b) ? (a) : (b))

void pdsptri( ivector, irange, n, dd, ee, dplus, lplus, lb, ub, ilb, iub, abstol,
	      meigval, vecZ, mapZ, eval, iscratch, iscsize,
	      dblptr, ibuffsize, scratch, ssize, info)
     Integer  *ivector, *irange, *n, *ilb, *iub, *meigval, *mapZ, *iscratch, *iscsize, *ibuffsize, *ssize, *info;
     DoublePrecision   *dd, *ee, *dplus, *lplus, *lb, *ub, *abstol, *eval, *scratch;
     DoublePrecision  **vecZ, **dblptr;
{
  
  /*
 *
 *  Purpose
 *  =======
 *
 *  pdsptri computes some or all of the eigenvalues and, optionally,
 *  eigenvectors of a real tri-diagonal eigenproblem of the form
 *
 *  A*x=(lambda) * x.
 *
 *  Here A is assumed to be symmetric, tri-diagonal.
 *
 *  Arguments
 *  =========
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
 *    nvecsZ = number of entries in mapZ equal to me
 *             (= count_list( me, mapZ, n ))
 *
 *-----------------------------------------------------------------
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
 *
 *  n       (input) INTEGER
 *          The number of rows and columns of the matrix A.
 *          N >= 0.
 *
 *  dd      (input) DOUBLE PRECISION array, dimension ( n )
 *          The diagonal elements of A.
 *
 *  ee      (input) DOUBLE PRECISION array, dimension ( n )
 *          The off-diagonal elements of A.  ee[0] is junk, the
 *          actual off-diagonal starts at ee[1].
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
 *          or equal to zero, then  EPS*|A|  will be used in its place,
 *          where |A| is the 1-norm of A.
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
 *          On entry, vecZ[i], i = 0 to nvecsZ-1, should point to a
 *          DOUBLE PRECISION array of length n.
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
 *                    i-th eigenvector, i = 0 to n-1 (current code
 *                    requires this even when computing only some eigenvalues.
 *          On exit:
 *            mapZ(i) = the id of the processor which actually owns the i-th
 *                      eigenvector, i = 0 to n-1.
 *
 *  eval    (output) DOUBLE PRECISION array, dimension (N)
 *          If INFO = 0, the eigenvalues  of the matrix
 *          in no particular order.
 *
 *  iscratch (workspace) INTEGER array,
 *           Length must be >= that returned from utility routine memreq_
 *
 *  iscsize  (input) INTEGER
 *           The number of usable elements in array "iscratch".
 *           Must be >= that returned from utility routine memreq_
 *
 *  dblptr   (workspace) DoublePrecision pointer to DoublePrecision (DoublePrecision **),
 *           Length must be >= that returned from utility routine memreq_
 *
 *  ibuffsize(input) INTEGER
 *           The number of usable elements in array "dblptr".
 *           Must be >= that returned from utility routine memreq_
 *
 *  scratch  (workspace) DOUBLE PRECISION array
 *           Length must be >= that returned from utility routine memreq_
 *
 *  ssize    (input) INTEGER
 *           The number of usable elements in array "scratch".
 *           Must be >= that returned in utility routine memreq_
 *
 *  INFO    (output) INTEGER
 *          = 0:  successful exit.
 *
 *          -50 <= INFO < 0, the -INFOth argument had an illegal value.
 *
 *          INFO = -51,      Input data which must be the same on all
 *                           processors is NOT the same on all processors
 *                           in mapZ.
 *
 *          = 1,             Error computing eigenvalues. PSTEBZ returned
 *                           a non-zero info, whose value was printed to
 *                           stderr.
 *
 *          = 2,             Error computing eigenvectors.  PSTEIN returned
 *                           a non-zero info, whose value was printed to
 *                           stderr.
 *
 *          Any processor with a negative INFO stops program execution.
 *
 *          If -50 <= INFO < 0, then bad data was passed to this routine
 *                              and no attempt is made to make sure that
 *                              INFO is the same on all processors.
 *        
 *          All other INFO      INFO should be the same on all processors in
 *                              mapZ.  If this routine calls a routine
 *                              which returns a negative info, then info
 *                              may not be the same on all processors.
 *                              This, however, should never happen.
 *---------------------------------------------------------------------
 */

/*
 * Local variables
 * ---------------
 */

    Integer             k, me, nn_proc, indx, msize, nsplit, neigval,
                        nproc, itmp, isize, nZ2, mapZ_0, nvecsZ, nvecsZ2,
                        linfo, maxinfo, i, j, num_procs;

    Integer            *i_scrat, *iblock, *isplit, *proclist;

    char                msg[ 35 ];
    char                msg2[ 35 ];

    Integer           **iptr;
    DoublePrecision   *d_scrat;
    DoublePrecision   **dd_scrat;
    
#ifdef TIMING
    extern DoublePrecision mxclock_();
    extern TIMINGG test_timing;
    DoublePrecision  t1, t2;
#endif

/*
 * External procedures
 * -------------------
 */

    extern Integer  count_list(), reduce_list2();
    extern void     xstop_(), pdiff(), pgexit();
    extern void     mdiff1_(), bbcast00();
    extern void     memreq_();

    extern void     pstebz_(), pstein();
    extern void     dshellsort_(), sorteig();

    extern Integer  mxmynd_(), mxnprc_();
    extern void     mxinit_();

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

   strcpy( msg,  "Error in pdsptri." );

#ifdef DEBUG1
   fprintf(stderr, "me = %d In pdsptri \n", me );
#endif

   /*
    * Set workspace pointers.
    */

   i_scrat  = iscratch;
   d_scrat  = scratch;
   dd_scrat = dblptr;

   /*
    *     Test the input parameters.
    */

   linfo = 0;

   if ( ivector == NULL )
      linfo = -1;
   else if ( irange == NULL )
      linfo = -2;
   else if ( n == NULL )
     linfo = -3;
   else if ( dd == NULL )
     linfo = -4;
   else if ( ee == NULL )
     linfo = -5;
   else if ( lb == NULL )
     linfo = -6;
   else if ( ub == NULL )
     linfo = -7;
   else if ( ilb == NULL )
      linfo = -8;
   else if ( iub == NULL )
      linfo = -9;
   else if ( abstol == NULL )
      linfo = -10;
   else if ( meigval == NULL )
      linfo = -11;
   else if ( vecZ == NULL )
      linfo = -12;
   else if ( mapZ == NULL )
      linfo = -13;
   else if ( eval == NULL )
      linfo = -14;
   else if ( iscratch == NULL )
      linfo = -15;
   else if ( iscsize == NULL )
      linfo = -16;
   else if ( dblptr == NULL )
      linfo = -17;
   else if ( ibuffsize == NULL )
      linfo = -18;
   else if ( scratch == NULL )
      linfo = -19;
   else if ( ssize == NULL )
      linfo = -20;
   else if ( info == NULL )
      linfo = -21;

   if ( linfo != 0 ) {
     if ( info != NULL )
        *info = linfo;

     fprintf( stderr, " %s me = %d argument %d is a pointer to NULL. \n",
              msg, me, -linfo );
     xstop_( &linfo );
     return;
   }

   *info    = 0;
   *meigval = 0;

   msize = *n;

   /*
    *  Quick return if possible.
    */

   if ( msize == 0 )
      return;

   /*
    *  Continue error checking.
    */

    itmp = 2;
    memreq_( &itmp, n, mapZ, mapZ, mapZ, &i, &j, &k, i_scrat );

   if ( *ivector < 0  || *ivector > 1 )
      *info = -1;
   else if ( *irange < 1  || *irange > 3 )
      *info = -2;
   else if ( *n < 0 )
      *info = -3;
   else if ( *irange == 2  && *lb >= *ub )
      *info = -7;
   else if ( *irange == 3  && *ilb < 1 )
      *info = -8;
   else if ( *irange == 3  && ( *iub < *ilb  ||  *iub > *n ) )
      *info = -9;
   else if ( *iscsize < i )
      *info = -16;
   else if ( *ibuffsize < k )
      *info = -18;
   else if ( *ssize < j )
      *info = -20;
   else

      for ( k = 0; k < *n; k++ )
         if ( mapZ[ k ] < 0  ||  mapZ[ k ] > nproc - 1 )
          *info = -13;
    
   if ( *info == 0 ) {

	nvecsZ = count_list( me, mapZ, &msize );

	if ( nvecsZ <= 0 )
	  return;

        for ( k = 0; k < nvecsZ; k++ )
	  if ( vecZ[ k ] == NULL )
	    *info = -12;
   }
    
   if ( *info != 0 ) {
      linfo = *info;
      fprintf( stderr, " %s me = %d argument %d has an illegal value. \n",
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
     *  Reduce mapZ to a single sorted list of processors.
     */
    
    proclist = i_scrat;
    nn_proc = reduce_list2( *n, mapZ, proclist);
    
    i_scrat += nn_proc;

    /*
     *  Check scaler Integer inputs.
     */
    
    i_scrat[ 0 ] = *ivector;
    i_scrat[ 1 ] = *irange;
    i_scrat[ 2 ] = *n;
    i_scrat[ 3 ] = *ilb;
    i_scrat[ 4 ] = *iub;
    
    isize = 5 * sizeof( Integer );
    strcpy( msg2,  "ivector,irange,n,ilb,or iub ");
    pdiff( &isize, (char *)i_scrat, proclist, &nn_proc, &i_scrat[5], msg, msg2, &linfo );
    
    pgexit( &linfo, msg, proclist, &nn_proc, d_scrat );
    
    if ( linfo != 0 ) {
	*info = -51;
	return;
    }

    /*
     *  Check other inputs.
     */
    
    maxinfo = 0;

    isize   = msize * sizeof( Integer );
    strcpy( msg2,  "mapZ ");
    pdiff( &isize, (char *)mapZ, proclist, &nn_proc, i_scrat, msg, msg2, &linfo );
    maxinfo = max( maxinfo, linfo );
    
    *d_scrat = *lb;
    *(d_scrat + 1) = *ub;
    *(d_scrat + 2) = *abstol;
    
    isize = 3 * sizeof( DoublePrecision );
    strcpy( msg2,  "lb,ub,or abstol ");
    pdiff( &isize, (char *)d_scrat, proclist, &nn_proc, i_scrat, msg, msg2, &linfo );
    maxinfo = max( maxinfo, linfo );
    
    linfo = maxinfo;
    
   pgexit( &linfo, msg, proclist, &nn_proc, d_scrat );

   if ( linfo != 0 ) {
      *info = -51;
      return;
   }

   /* ----------------------------------------
    * All input data is good. Start computing. 
    * ----------------------------------------
    */

    /*
     * Release proclist.
     */

    proclist = NULL;
    i_scrat -= nn_proc;

    /*
     * Initialize workspace.
     */

    iblock = i_scrat;
    i_scrat += msize;

    isplit = i_scrat;
    i_scrat += msize;


    iptr = (Integer **) dd_scrat;
    

    mapZ_0 = *mapZ;
    
    /*
     * Compute eigenvalues.
     */
    
#ifdef TIMING
	t1 = mxclock_();
#endif
	
	pstebz_( irange, &msize, lb, ub, ilb, iub, abstol, dd, ee,
		 dplus, lplus,
		 mapZ, &neigval, &nsplit, eval, iblock, isplit,
		 d_scrat, i_scrat, info);
	
#ifdef TIMING
	mxsync_();
	t2 = mxclock_();
	test_timing.pstebz = t2 - t1;
#endif

    *meigval = neigval;

    if( *info != 0 ) {
      fprintf(stderr, " %s me = %d pstebz returned info = %d \n", msg, me, *info );
      *info = 1;
      return;
    }
    
    /*
     * Sort eigenvalues and return if there aren't any eigenvectors
     * to be computed.
     */
    
    if ( neigval == 0 )
      return;
    
    if ( *ivector == 0 ) {
      dshellsort_( &neigval, eval );
      return;
    }
    
    /*
     * Compute eigenvectors.
     */
    
    nvecsZ2 = count_list( me, mapZ, meigval );
    if (nvecsZ2 > 0) {
      
#ifdef TIMING
      t1 = mxclock_();
#endif
      
      pstein( &msize, dd, ee, &neigval, eval, iblock, &nsplit, isplit,
	     mapZ, vecZ, d_scrat, i_scrat, iptr, info);
      
#ifdef TIMING
      mxsync_();
      t2 = mxclock_();
      test_timing.pstein = t2 - t1;
#endif
      
      if( *info != 0 ) {
        fprintf(stderr, " %s me = %d pstein returned info = %d \n", msg, me, *info );
        *info = 2;
      }
    }
    
    /*
     * Release iblock and isplit.
     */

    iblock = NULL;
    i_scrat -= msize;

    isplit = NULL;
    i_scrat -= msize;
    
    /*
     * Send info, mapZ to any processors in 
     * {mapZ[neigval:*n-1]}) - {mapZ[0:neigval-1]}
     */
    
    nZ2 = *n - neigval;
    mdiff1_( &nZ2, &mapZ[neigval], meigval, mapZ, i_scrat, &num_procs ); 
    
    if( num_procs > 0 ){ /* make sure everyone has the same maps */
      i_scrat[num_procs] = mapZ_0;
      num_procs++;
      bbcast00( (char *)info, sizeof(Integer),       2, mapZ_0, num_procs, i_scrat );
      bbcast00( (char *)mapZ, msize*sizeof(Integer), 3, mapZ_0, num_procs, i_scrat);
    }
    
    for ( indx = neigval; indx < msize; indx++ )
      mapZ[indx] = -1;
    
    if ( *info != 0 )
      return;
    
    /*
     * Sort eigenvalues, mapZ and eigenvectors
     */
    
    sorteig( &msize, &neigval, vecZ, mapZ, eval, iscratch, scratch );
    
#ifdef DEBUG1
    fprintf(stderr, "me = %d Exiting pdsptri \n", me );
#endif
     
    return;
  }
