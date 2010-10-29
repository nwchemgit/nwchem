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

#include "globalp.c.h"
#include "blas_lapack.h"


#define min(a,b) ((a) < (b) ? (a) : (b))
#define max(a,b) ((a) > (b) ? (a) : (b))

#define ZERO  ((DoublePrecision) 0.0e0)

extern DoublePrecision psigma, psgn;

void pstebz9_( job, n, lb, ub, jjjlb, jjjub, abstol, d, e,
		dplus, lplus, mapZ, neigval,
		nsplit, eval, iblock, isplit, work, iwork, info)
     Integer            *job, *n, *jjjlb, *jjjub, *mapZ, *neigval, *nsplit,
  *iblock, *isplit, *iwork, *info;
     DoublePrecision         *lb, *ub, *abstol, *d, *e, *dplus, *lplus, *eval, *work;
{
  
/*
 *  Parallel version of LAPACK's DSTEBZ.
 *
 *  This routine is directly callable from both FORTRAN and C.
 *  The documentation below always uses FORTRAN array indexing,
 *  i.e., 1 to N, rather then C array indexing, i.e.,  0 to N-1.
 *  This should be kept in mind when calling this routine from C.
 *  Also when calling this routine from C 
 *    INTEGER (array)           means "pointer to Integer" and
 *    DOUBLE PREICISION (array) means "pointer to DoublePrecision"
 *
 *
 *  Purpose
 *  =======
 *
 *  PSTEBZ computes the eigenvalues of a symmetric tridiagonal
 *  matrix T.  The user may ask for all eigenvalues, all eigenvalues
 *  in the half-open interval (LB, UB], or the JJJLB-th through JJJUB-th
 *  eigenvalues.
 *
 *  To avoid overflow, the matrix must be scaled so that its
 *  largest element is no greater than overflow**(1/2) *
 *  underflow**(1/4) in absolute value, and for greatest
 *  accuracy, it should not be much smaller than that.
 *
 *  See W. Kahan "Accurate Eigenvalues of a Symmetric Tridiagonal
 *  Matrix", Report CS41, Computer Science Dept., Stanford
 *  University, July 21, 1966.
 *
 *  Arguments
 *  =========
 *
 *
 *  JOB     (input) INTEGER
 *          = 1  : ("All")   all eigenvalues will be found.
 *          = 2  : ("Value") all eigenvalues in the half-open interval
 *                           (LB, UB] will be found.
 *          = 3  : ("Index") the JJJLB-th through JJJUB-th eigenvalues (of the
 *                           entire matrix) will be found.
 *
 *  N       (input) INTEGER
 *          The order of the tridiagonal matrix T.  N >= 0.
 *
 *  LB      (input) DOUBLE PRECISION
 *  UB      (input) DOUBLE PRECISION
 *          If JOB=2, the lower and upper bounds of the interval to
 *          be searched for eigenvalues.  Eigenvalues less than or equal
 *          to LB, or greater than UB, will not be returned.  LB < UB.
 *          Not referenced if JOB = 1 or 3.
 *
 *  JJJLB   (input) INTEGER
 *  JJJUB   (input) INTEGER
 *          If JOB=3, the indices (in ascending order) of the
 *          smallest and largest eigenvalues to be returned.
 *          1 <= JJJLB <= JJJUB <= N, if N > 0.
 *          Not referenced if JOB = 1 or 2.
 *
 *  ABSTOL  (input) DOUBLE PRECISION
 *          The absolute tolerance for the eigenvalues.  An eigenvalue
 *          (or cluster) is considered to be located if it has been
 *          determined to lie in an interval whose width is ABSTOL or
 *          less.  If ABSTOL is less than or equal to zero, then ULP*|T|
 *          will be used, where |T| means the 1-norm of T.
 *
 *          Eigenvalues will be computed most accurately when ABSTOL is
 *          set to twice the underflow threshold 2*DLAMCH('S'), not zero.
 *
 *  D       (input) DOUBLE PRECISION array, dimension (N)
 *          The n diagonal elements of the tridiagonal matrix T.
 *
 *  E       (input) DOUBLE PRECISION array, dimension (N)
 *          The first element of E, E(1), is junk, the rest of E, E(2:N),
 *          contains the (n-1) off-diagonal elements of the tridiagonal
 *          matrix T.
 *
 *  MAPZ    (input) INTEGER array, dimension (N)
 *          A list of the ids of the processors which are to participate
 *          in the eigenvalue computation.
 *
 *  NEIGVAL (output) INTEGER
 *          The actual number of eigenvalues found. 0 <= NEIGVAL <= N.
 *
 *  NSPLIT  (output) INTEGER
 *          The number of diagonal blocks in the matrix T.
 *          1 <= NSPLIT <= N.
 *
 *  EVAL    (output) DOUBLE PRECISION array, dimension (N)
 *          On exit, the first NEIGVAL elements of EVAL will contain the
 *          eigenvalues.  (PSTEBZ may use the remaining N-NEIGVAL elements as
 *          workspace.)
 *          The eigenvalues will be grouped by split-off block (see IBLOCK,
 *          ISPLIT) and ordered from smallest to largest within the block.
 *
 *  IBLOCK  (output) INTEGER array, dimension (N)
 *          At each row/column j where E(j) is zero or small, the
 *          matrix T is considered to split into a block diagonal
 *          matrix.  On exit, if INFO = 0, IBLOCK(i) specifies to which
 *          block (from 1 to the number of blocks) the eigenvalue EVAL(i)
 *          belongs.  (DSTEBZ may use the remaining N-NEIGVAL elements as
 *          workspace.)
 *
 *  ISPLIT  (output) INTEGER array, dimension (N)
 *          The splitting points, at which T breaks up into submatrices.
 *          The first submatrix consists of rows/columns 1 to ISPLIT(1),
 *          the second of rows/columns ISPLIT(1)+1 through ISPLIT(2),
 *          etc., and the NSPLIT-th consists of rows/columns
 *          ISPLIT(NSPLIT-1)+1 through ISPLIT(NSPLIT)=N.
 *          (Only the first NSPLIT elements will actually be used, but
 *          since the user cannot know a priori what value NSPLIT will
 *          have, N words must be reserved for ISPLIT.)
 *
 *  WORK    (workspace) DOUBLE PRECISION array, dimension ( )
 *
 *  IWORK   (workspace) INTEGER array, dimension ( )
 *
 *  INFO    (output) INTEGER
 *
 *          PSTEBZ attempts to return the same INFO on all processors in MAPZ.
 *          Currently, however, if the input data is invalid, -50 < INFO < 0,
 *          then INFO will be different on different processors.
 *
 *          = 0:   successful exit
 *
 *          < 0 &
 *          > -50: if INFO = -i, the i-th argument had an illegal value
 *
 *          = -51: Then the input data was not the same on all processors in MAPZ
 *
 *          > 0:   some or all of the eigenvalues failed to converge or
 *                 were not computed, or they were computed incorrectly:
 *
 *                 =1: routine DSTEBZ3 returned a non-zero info.  Meaning
 *                     that on one or more processors some or all of the
 *                     eigenvalues failed to converge or were not computed
 *                     In this case the processor with a non-zero info from
 *                     DSETBZ3 prints and error message to stderr.
 *
 *                 =2: NSPLIT and/or ISPLIT were not the same on all processors
 *
 *                 =3: Relative processor i computed an eigenvalue larger than
 *                     the smallest eigenvalue computed by relative processor i+1
 *                     for some i.  This should not occur using the current algorithms.
 *
 *                 =4: The number of eigenvalues found in a block of T is bigger
 *                     then the dimension of the block.
 *
 *             
 *               In theory INFO = 2, 3, or 4 should never occur.  If they do occur
 *               then DSTEBZ3 failed to correctly compute the requested eigenvalues.
 *
 *               Note that DSTEBZ3 is a modification of LAPACK's DSTEBZ.  The
 *               modifications had to be made to avoid gettings INFOs like 3 and 4.
 *
 */

   /*
    *  Local Variables
    */

   static Integer      INT = 10, INT2 = 20, DOUBLE = 200;

   char msg[35], *cptr;
   char msg2[35];
   Integer range, order;

   Integer             indx, il, iu, ifakeme, msgli, msglr, itype,
                   nhigh, numeig, irem, isize, ival, linfo, isize1,
                   iii, m, nlow, me, ncol, ii, junk, k, nn_procs,
                   nproc, msize, maxinfo, *iptr, *i_work, *proclist, ncols;

   DoublePrecision         lstmax, emax, emin, leig, reig, eps;
   Integer i, i1split, jsplit, blksz, j, jjj;
   DoublePrecision *dptr, *lptr, onenrm, dummy;
   
   /*
    *  External Procedures
    */

   extern Integer      mxmynd_(), mxnprc_(), mxwrit_(), mxread_(), mxbrod_();

   extern void     sort_();
   extern void     dstebz3_();

   extern Integer  menode_();
   extern Integer  neblw2_();
   extern Integer  mapchk_();
   extern void     xstop_(), pdiff(), pgexit();

   
/*
 *  ---------------------------------------------------------------
 *                      Executable Statements
 *  ---------------------------------------------------------------
 */
   
   /*
    *  Get this processor's id, and the number of allocated nodes.
    */
   
   me    = mxmynd_();
   nproc = mxnprc_();

   strcpy( msg,  "Error in pstebz." );
   
#ifdef DEBUG1
   fprintf(stderr, "me = %d In pstebz \n", me );
#endif
   
   /*
    *     Test the input parameters.
    */
   
   linfo = 0;

   if ( job == NULL )
      linfo = -1;
   else if ( n == NULL )
      linfo = -2;
   else if ( lb == NULL )
      linfo = -3;
   else if ( ub == NULL )
      linfo = -4;
   else if ( jjjlb == NULL )
      linfo = -5;
   else if ( jjjub == NULL )
      linfo = -6;
   else if ( abstol == NULL )
      linfo = -7;
   else if ( d == NULL )
      linfo = -8;
   else if ( e == NULL )
      linfo = -9;
   else if ( mapZ == NULL )
      linfo = -10;
   else if ( neigval == NULL )
      linfo = -11;
   else if ( nsplit == NULL )
      linfo = -12;
   else if ( eval == NULL )
      linfo = -13;
   else if ( iblock == NULL )
      linfo = -14;
   else if ( isplit == NULL )
      linfo = -15;
   else if ( work == NULL )
      linfo = -16;
   else if ( iwork == NULL )
      linfo = -17;
   else if ( info == NULL )
      linfo = -18;

   if ( linfo != 0 ) {
     if ( info != NULL )
        *info = linfo;

     fprintf( stderr, " %s me = %d argument %d is a pointer to NULL. \n",
              msg, me, -linfo );
     xstop_( &linfo );
     return;
   }
   
   *info = 0;
   *neigval = 0;
   *nsplit  = 0;

    msize = *n;

   /*
    *  Quick Return if possible.
    */

   if ( *n == 0 )
      return;

   /*
    *  Continue error checking.
    */

   if ( *job < 1  || *job > 3 )
     linfo = -1;
   else if ( *n < 0 )
      linfo = -2;
   else if ( *job == 2  && *lb >= *ub )
      linfo = -4;
   else if ( *job == 3  && *jjjlb < 1 )
      linfo = -5;
   else if ( *job == 3  && ( *jjjub < *jjjlb  ||  *jjjub > *n ) )
      linfo = -6;
   else if( mapchk_( n, mapZ ) != 0 )
      linfo = -10;
   
   if ( linfo != 0 ) {
       fprintf( stderr, " %s me = %d argument %d has an illegal value. \n",
                msg, me, -linfo );
       *info = linfo;
       xstop_( info );
       return;
   }
   
   /*
    *  ------------------------------------------------
    *  No local errors, compare data across processors.
    *  ------------------------------------------------
    */
   
   /*
    *  Determine the number of unique processor ids in mapZ, nn_procs,
    *  and this processors relative position in mapZ.
    */
       
    i_work = iwork;

    proclist = i_work;
    nn_procs = reduce_list2( *n, mapZ, proclist );
    i_work += nn_procs;

    ifakeme = menode_(&nn_procs, proclist );
       
    if (ifakeme < 0)
      return;
       
   /*
    *  Check Integer scaler inputs, mainly interested in *n.
    *
    *  If there are any difference then exit since n may be different
    *  in which case there is no point in checking arrays of length n.
    */
   
   i_work[ 0 ] = *job;
   i_work[ 1 ] = *n;
   i_work[ 2 ] = *jjjlb;
   i_work[ 3 ] = *jjjub;
   
   isize = 4 * sizeof(Integer);
   strcpy( msg2,  "job,n,jjlb,or jjub " );
   pdiff( &isize, (char *) i_work, proclist, &nn_procs, i_work+4, msg, msg2, &linfo );

   pgexit( &linfo, msg, proclist, &nn_procs, work );

   if ( linfo != 0 ) {
      *info = -51;
      return;
   }

   /*
    *  Check remaining inputs.
    */

   maxinfo = 0;

   isize = msize * sizeof(Integer);
   strcpy( msg2,  "mapZ " );
   pdiff( &isize, (char *) mapZ, proclist, &nn_procs, i_work, msg, msg2, &linfo );
   maxinfo = max( maxinfo, linfo );

   work[ 0 ] = *lb;
   work[ 1 ] = *ub;
   work[ 2 ] = *abstol;
   
   isize = 3 * sizeof(DoublePrecision);
   strcpy( msg2,  "lb,ub,abstol " );
   pdiff( &isize, (char *) work, proclist, &nn_procs, (Integer *) (work+3), msg, msg2, &linfo );
   maxinfo = max( maxinfo, linfo );

   isize = msize * sizeof(DoublePrecision);
   strcpy( msg2,  "d " );
   pdiff( &isize, (char *) d, proclist, &nn_procs, (Integer *) work, msg, msg2, &linfo );
   maxinfo = max( maxinfo, linfo );

   strcpy( msg2,  "e " );
   pdiff( &isize, (char *) e, proclist, &nn_procs, (Integer *) work, msg, msg2, &linfo );
   maxinfo = max( maxinfo, linfo );

   pgexit( &linfo, msg, proclist, &nn_procs, work );

   if ( linfo != 0 ) {
      *info = -51;
      return;
   }

   /* ----------------------------------------
    * All input data is good. Start computing. 
    * ----------------------------------------
    */

   /*
    *  Compute the index of the first and last eigenvalues to be found.
    */

   if ( *job == 1) {
      nlow  = 1;
      nhigh = *n;
   }
   else if ( *job == 3 ) {
      nlow  = *jjjlb;
      nhigh = *jjjub;
   }
   else {

      /*
       *  Must compute nlow and nhigh.  To guarantee that all processors
       *  get the same values we compute these numbers on processor
       *  ifakeme = 0 and broadcast the results to everyone else.
       */

      if ( ifakeme == 0) {
         nlow  = 1 + neblw2_( n, lb, d, e+1, work, &linfo );
         nhigh =     neblw2_( n, ub, d, e+1, work, &linfo );

         i_work[ 0 ] = nlow;
         i_work[ 1 ] = nhigh;
      }

      if ( nn_procs > 1) {
         isize = 2 * sizeof( Integer );
         itype = 15;
         mxbrod_( i_work, proclist, &isize, &nn_procs, proclist, &itype );

         nlow  = i_work[ 0 ];
         nhigh = i_work[ 1 ];
      }
   }

   /*
    * Have each processor compute "ncol" eigenvalues, though 
    * "irem" processors have to compute "ncol+1" eigenvalues.
    */

   numeig   = nhigh - nlow + 1;
   *neigval = numeig;
   
   if ( numeig == 0 )
     return;
   
   ncol = numeig / nn_procs;
   irem = (numeig % nn_procs);

/*    printf("got here in pdspevx 9 "); */

   goto JUNK;
   
   if ( ifakeme < irem ) {
       il = ifakeme * (ncol + 1) + nlow;
       iu = il + ncol;
   }
   else {
       il = irem + ifakeme * ncol + nlow;
       iu = il + ncol - 1;
   }
   
   if ( iu >= il ) {
     /*
     range = 'I';
     order = 'B';
     */
     range = 3;
     order = 1;
     m = 0;
     dstebz3_( &range, &order, n, lb, ub, &il, &iu, abstol, d, e+1,
	       &m, nsplit, eval, iblock, isplit, work,
	       i_work, info);
     
     if ( *info != 0 ) {
       fprintf(stderr, " me = %d ifakeme=%d dstebz3 returned info = %d to pstebz \n",
	       me, ifakeme, *info );
       fprintf(stderr, " me = %d ifakeme=%d dstebz3 also returned neig= %d il=%d iu=%d \n",
	       me, ifakeme, m, il, iu );
       
       *info = 1;
     }
     
   }
   
   /*
    *  Broadcast nsplit and isplit, and make sure they are the same
    *  on all processors which computed eigenvalues.  This is just paranoia.
    */
   
   if ( ifakeme == 0  ||  iu < il ) {
     isize = sizeof( Integer );
     itype = 16;
     mxbrod_( nsplit, proclist, &isize, &nn_procs, proclist, &itype );

     isize = *nsplit * sizeof( Integer );
     itype = 17;
     mxbrod_( isplit, proclist, &isize, &nn_procs, proclist, &itype );
   }
   else {

     isize = sizeof( Integer );
     itype = 16;
     mxbrod_( i_work, proclist, &isize, &nn_procs, proclist, &itype );
     
     if ( *info == 0  &&  *nsplit != i_work[ 0 ]  ) {
       fprintf(stderr, " me = %d in pstebz. My nsplit= %d differs from others= %d \n",
               me, *nsplit, i_work[0] );

       *info = 2;
       *nsplit = i_work[ 0 ];
     }

     isize = *nsplit * sizeof( Integer );
     itype = 17;
     mxbrod_( i_work, proclist, &isize, &nn_procs, proclist, &itype );

     if ( *info== 0  ) {

       linfo = 0;
       for (indx = 0; indx < *nsplit; indx++)
         if ( isplit[ indx ] != i_work[ indx ] )
           linfo = 1;

       if ( linfo != 0 ) {
         fprintf(stderr, " me = %d in pstebz. My isplit differs from others \n", me );

         *info = 2;
       }
     }
   }
   
   /*
    *  Send info, iblock and eval to processor ifakeme = 0.
    */
   
   isize1  = sizeof(Integer);

   itype = 18;

   if ( ifakeme > 0  &&  iu >= il ) {
       cptr = (char *) ( &junk );
       isize = mxread_( cptr, &isize1, &proclist[0], &itype );
       
       cptr = (char *) ( info );
       isize = mxwrit_( cptr, &isize1, &proclist[0], &INT2 );
       
       msgli = (iu - il + 1) * sizeof(Integer);
       msglr = (iu - il + 1) * sizeof(DoublePrecision);
       
       cptr = (char *) iblock;
       isize = mxwrit_( cptr, &msgli,  &proclist[0], &INT );
       
       cptr = (char *) eval;
       isize = mxwrit_( cptr,   &msglr,  &proclist[0], &DOUBLE );
   }
   
   /*
    * Relative processor 0 receives
    */

   if ( ifakeme == 0  && nn_procs > 1 ) {
      junk = 0;

      /* Find processor ifakeme = 0's maximum eigenvalue. */

      lstmax = eval[0];
      for (indx = 1; indx < iu-il+1; indx++)
        lstmax = max( lstmax, eval[indx] );

      ncols = iu-il+1;
      iii   = ncols;
      for (indx = 1; indx < nn_procs; indx++) {

         if ( indx == irem ) {
           ncols = ncol;
           if( ncols == 0 )
             break;
         }

	 cptr = (char *) &junk;
	 isize = mxwrit_( cptr, &isize1, &proclist[ indx ], &itype );
	   
	 cptr = (char *) &linfo;
	 isize = mxread_( cptr, &isize1, &proclist[ indx ], &INT2 );

	 if( linfo != 0 )  {
           if( *info == 0 )
	     *info = linfo;
           else
	     *info = min( *info, linfo );
         }
	   
	 msgli = ncols * sizeof(Integer);
	 msglr = ncols * sizeof(DoublePrecision);
	   
	 cptr = (char *) ( iblock + iii );
	 isize = mxread_( cptr, &msgli, &proclist[ indx ], &INT );
	   
	 cptr = (char *) ( eval + iii );
	 isize = mxread_( cptr,   &msglr, &proclist[ indx ], &DOUBLE );

         /* make sure minimum eigenvalue computed by current processor
          * is >= the largest eigenvalue computed by previous processors.
          */

         emax = eval[iii];
         emin = eval[iii];
         for (ii = iii+1; ii < iii+ncols; ii++) {
           emax = max( emax, eval[ii] );
           emin = min( emin, eval[ii] );
         }

         if( emin < lstmax && *info == 0 ) {
           fprintf( stderr, " ERROR in pstebz, computed eigenvalues are not monotonic.\n");
           *info = 3;
         }

         lstmax = emax;

         iii += ncols;
      }
   }

   /*
    * Make sure we didn't find too may eigenvalues in any of
    * the blocks.
    */

   if ( *info != 0 ) {
/*
	if parallel multisection gives errors then do it serially
*/
     range = 3;
     order = 1;
     m = 0;
     il = 1;
     iu = *n;
     dstebz3_( &range, &order, n, lb, ub, &il, &iu, abstol, d, e+1,
	       &m, nsplit, eval, iblock, isplit, work,
	       i_work, info);
    }
    
    
    if ( ifakeme == 0  &&  *info == 0 ) {

      for ( indx = 0; indx < *nsplit; indx++ )
         i_work[ indx ] = 0;

      for ( indx = 0; indx < *neigval; indx++ )
         (i_work[ iblock[ indx ] - 1 ])++;

      if ( i_work[ 0 ] > isplit[ 0 ] ) {
         *info = 4;
         fprintf( stderr, " ERROR in pstebz, too many eigenvalues in block %d\n", 0 );
      }

      for ( indx = 1; indx < *nsplit; indx++ )
         if ( i_work[ indx ] > isplit[ indx ] - isplit[ indx - 1 ] ) {
           fprintf( stderr, " ERROR in pstebz, too many eigenvalues in block %d\n", indx );
           *info = 4;
         }

      if( nn_procs > 1  &&  *info == 0 )
	sort_( n, &numeig, nsplit, isplit, iblock, i_work, eval, work );
    }
    
    /*
     *  Broadcast info, iblock, and eval to all processors in mapZ.
     */
    
    mxbrod_( info, proclist, &isize1, &nn_procs, proclist, &INT);
 
    itype = 19;
    isize = numeig * sizeof(Integer);
    mxbrod_( iblock, proclist, &isize, &nn_procs, proclist, &itype);
    
    isize = numeig * sizeof(DoublePrecision);
    mxbrod_( (Integer *) eval, proclist, &isize, &nn_procs, proclist, &DOUBLE);

JUNK:
    
    range = 3;
    order = 1;
    m = 0;
    il = 1;
    iu = *n;
    dstebz3_( &range, &order, n, lb, ub, &il, &iu, abstol, d, &e[1],
	      &m, nsplit, eval, iblock, isplit, work,
	      i_work, info);
    
    sort_( n, &numeig, nsplit, isplit, iblock, i_work, eval, work );
    
    leig = eval[0]; /* left most e-val */
    psgn = 1.0;
    psigma = 0.;
    eps = DLAMCHE;
    
    
    /*
      for ( indx = 0; indx < *n; indx++ )
      work[indx] = d[indx];
      
      for ( indx = 0; indx < *n; indx++ )
      work[msize + indx] = e[indx];
      
      dsterf_(&msize, work, &work[msize], info);
      
      for ( indx = 0; indx < *n; indx++ )
      eval[indx] = work[indx];
      
      *neigval = msize;
      leig = work[0];
      */
    
    if ( leig <= 0. ){
      psgn = 1.;
      psigma = -(fabs(leig)+sqrt(eps));
    }
    
    for ( indx = 0; indx < *neigval; indx++ ){
      eval[indx] -= psigma;
    }
    
    
    for (i = 0; i < msize; i++ )
      work[i] = d[i] - psgn*psigma;
    
    i1split = 0;
    for ( iii = 0; iii < *nsplit; iii++ ){
      jsplit = isplit[iii];
      blksz = jsplit-i1split;
      dptr = &dplus[i1split];
      lptr = &lplus[i1split];
      /*
	LDL' factorization of tridiagonal
	*/
      peigs_tldlfact(&blksz, &work[i1split], &e[i1split], dptr, lptr);
      i1split = jsplit;
    }
    
    
#ifdef DEBUG1
    printf(" me = %d ifakeme = %d Exiting pstebz \n", me, ifakeme);
    for ( iii = 0; iii < *nsplit; iii++ ){
      printf(" me = %d iii  %d iblock = %d isplit %d \n", me, iii, iblock[iii], isplit[iii]);
    }
#endif
    /*
      printf(" exiting pstebz9.c \n");
      */
    
   return;
}
