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
 *  -- PEIGS  routine (version 3.0) --
 *     Pacific Northwest Laboratory
 *     July 28, 1995
 *
 *======================================================================
 */
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "globalp.c.h"
#include "blas_lapack.h"

#define min(a,b) ((a) < (b) ? (a) : (b))
#define max(a,b) ((a) > (b) ? (a) : (b))
#define ffabs(a) ((a) > (0.) ? (a) : (-a))

#define ZERO  ((DoublePrecision) 0.0e0)

extern DoublePrecision psigma, psgn, peigs_scale, peigs_shift;

void pscale_( job, n, lb, ub, jjjlb, jjjub, abstol,
	     d, e, dplus, lplus, mapZ, neigval,
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

 bug fixed with new lapack code

 *
 */

   /*
    *  Local Variables
    */

  /*
   static Integer      INT = 10, INT2 = 20, DOUBLE = 200, IONE=1;   
  */

   char msg[35];

   /*
   char msg2[35];
   Integer range, order;
   */
   
   Integer         il, linfo, me,
     nproc, msize, maxinfo; 
   Integer i, idummy;
   DoublePrecision onenrm;
   
   /*
   DoublePrecision         lstmax, emax, emin, ulp, safemn, onenrm;
   DoublePrecision         peigs_leig, peigs_reig;
   */
   
   /*
    *  External Procedures
    */

   extern Integer      mxmynd_(), mxnprc_(), mxwrit_(), mxread_(), mxbrod_();

   extern void     sort_();
   extern void     dstebz_();

   extern Integer  menode_();
   extern Integer  neblw2_();
   extern Integer  mapchk_();
   extern void     xstop_(), pdiff(), pgexit();
   DoublePrecision dummy, tmp; 
   extern DoublePrecision dlamch_();
   

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

   *info = 0;
   *neigval = 0;
   *nsplit  = 0;

    msize = *n;
    
   /*
    *  Quick Return if possible.
    */

   /*
    *  Check remaining inputs.
    */

   maxinfo = 0;
   
   msize = *n;
   
   onenrm = ffabs( d[0] ) + ffabs( e[1] );
   for (i = 1; i < msize-1; i++) {
     tmp = ffabs(d[i]) + ffabs(e[i]) + ffabs(e[i + 1]);
     onenrm = max(onenrm, tmp);
   }
   tmp = ffabs(d[msize-1]) + ffabs(e[msize-1]);
   onenrm = max(onenrm, tmp);

   idummy = 0.0e0;
   peigs_scale =DLAMCHS;
   for (i = 0; i < msize; ++i) {
     peigs_scale = max(peigs_scale, ffabs(d[i]));
     if ( d[i] < 0.0e0 ){
       idummy = -1;
     }
   }

   for (i = 1; i < msize; ++i)
     peigs_scale = max(peigs_scale, ffabs(e[i]));
   
   peigs_shift = 0.0e0;
   if ( idummy == -1 ) {
     peigs_shift = -onenrm;
     peigs_scale = peigs_scale + onenrm;
   }
   
   /*
     il = 1;
     iu = 1;
     range = 3;
     order = 1;
     m = 0;
     *info = 0;
     
     dstebz3_( &range, &order, n, lb, ub, &il, &iu, abstol, d, e+1,
     &m, nsplit, eval, iblock, isplit, work, iwork, &linfo);
     if ( *info != 0 ) {
     printf(" error in stebz3 %d info %d leig %g  \n", me, *info, leig);
     }
     peigs_leig = eval[0];
     
     il = msize;
     iu = msize;
     range = 3;
     order = 1;
     m = 0;
     *info = 0;
     dstebz3_( &range, &order, n, lb, ub, &il, &iu, abstol, d, e+1,
     &m, nsplit, eval, iblock, isplit, work, iwork, &linfo);
     
     if ( *info != 0 ) {
     printf(" error in stebz3 %d info %d leig %g  \n", me, *info, reig);
     }
     
     peigs_reig = eval[0];
     peigs_leig -= DLAMCHE;
     peigs_reig += DLAMCHE;
     peigs_scale = peigs_reig - peigs_leig;
     eps = DLAMCHE;
   */
   
   /*
     if ( ffabs( peigs_leig ) <= ffabs( peigs_reig )) {
     peigs_shift = peigs_leig - 2.0e0*eps*onenrm + onenrm;
     peigs_scale = 1.0e0;
     }
     else
     peigs_shift = peigs_reig + 2.0e0 * eps * onenrm + onenrm;
     peigs_scale = -1.0e0;
     }
     printf(" peigs_leig %g peigs_reig %g \n", peigs_leig, peigs_reig );
   */
   

   /*
     printf(" peigs_shift %g peigs_scale %g \n", peigs_shift, peigs_scale );
     */
   
   dummy = 1.0e0/peigs_scale;   
   for ( il=0; il < msize; il++ ){
     d[il] = dummy*(d[il] - peigs_shift );
   }
   
   
   for ( il = 1; il < msize; il++ ){
     e[il] *= dummy;
   }
   
   *info = 0;
   
 }
   
