/*
 $Id$
 *======================================================================dstebz
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
#include <string.h>

#include "globalp.c.h"
#include "blas_lapack.h"

#define min(a,b) ((a) < (b) ? (a) : (b))
#define max(a,b) ((a) > (b) ? (a) : (b))
#define ffabs(a) ((a) > (0.) ? (a) : (-a))

#define ZERO  ((DoublePrecision) 0.0e0)

extern DoublePrecision psigma, psgn;

void pstebz10_( job, n, lb, ub, jjjlb, jjjub, abstol, d, e, dplus, lplus, mapZ, neigval,
		nsplit, eval, iblock, isplit, work, iwork, info)
     
     Integer            *job, *n, *jjjlb, *jjjub, *mapZ, *neigval, *nsplit,
  *iblock, *isplit, *iwork, *info;
     DoublePrecision         *lb, *ub, *abstol, *d, *e, *dplus, *lplus, *eval, *work;
{
  
  /*
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

  char msg[35];
  char msg2[35];
  Integer range, order;
  
  Integer         il, iu, ifakeme, itype,
    nhigh, numeig, isize, linfo,
    iii, m, nlow, me, nn_procs,
    nproc, msize, maxinfo, *i_work, *proclist;
  Integer j, i, i1split, jsplit, blksz;
  DoublePrecision *dptr, *lptr;
  
  DoublePrecision ulp, safemn, onenrm;
  
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
   DoublePrecision eps, tmp;
   extern void dsterf_();
   extern void peigs_tldlfact();
   extern DoublePrecision dlamch_();

   FILE *file;
   char filename[40];

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

   sprintf( filename, "pdspevx.%d", (int)me);
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
              msg, (int)me, (int)-linfo );
     xstop_( &linfo );
     return;
   }
   
   *info = 0;
   *neigval = 0;
   *nsplit  = 0;
   
    msize = *n;

   if ( msize == 1 ) {
     *nsplit = 1;
     isplit[0] = 1;
     if ( *lb > d[0] || *ub < d[0])
       *neigval = 0;
     else {
       eval[0] = d[0];
       iblock[0] = 1;
       *neigval = 1;
     }
     return;
   }
   
   onenrm = ffabs( d[0] ) + ffabs( e[1] );
   for (i = 1; i < msize-1; ++i) {
     tmp = ffabs(d[i]) + ffabs(e[i]) + ffabs(e[i + 1]);
     onenrm = max(onenrm, tmp);
   }
   tmp = ffabs(d[msize-1]) + ffabs(e[msize-1]);
   onenrm = max(onenrm, tmp);
   
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
                msg, (int)me, (int)-linfo );
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
   *neigval = *n;
   
   if ( numeig == 0 )
     return;
#ifdef DEBUG2
   for ( m = 0; m < *n; m++ )
     printf(" on entrance to pstebz10 i = %d d %g e %g \n", m, d[m], e[m]);
#endif


   
#ifdef DEBUG3
   for ( il = 0; il < m;  il++ )
     printf(" il = %d eval %f info = %d \n", il, eval[il], *info);
   
#endif
   

   
   
  /*
    computing split
    */
  
   eps = 2.0*sqrt(DLAMCHE);
   ulp = DLAMCHU;
   safemn = DLAMCHS;
   
   if ( *n == 1 ){
     *nsplit = 1;
     i1split = 1;
     isplit[0] = 1;
     eval[0] = d[0];
     m = 1;
     iblock[0]= 1;
     psgn = 1.0;
     psigma = 0.;
     return;
   }
   
   /*
     eps = 2.0*DLAMCHE;
     shift = 2.*eps*max(ffabs(leig), ffabs(reig));
     if ( ffabs(leig) <= ffabs(reig) ){
     psigma = leig -  shift;
     psgn = 1.0;
     }
     else {
     psigma = reig +  shift;
     psgn = -1.0;
     }
   */
   
   
   /*
     if ( ffabs(leig) <= ffabs(reig) ){
     psigma = reig -  shift;
     psgn = 1.0;
     }
     else {
     psigma = reig +  shift;
     psgn = -1.0;
     
     }
     */
   

   /*
     can be bad
   */
   /*
     leig = eval[0]; 
     shift = d[0] - fabs(e[1]);
     for ( i = 1; i < msize - 1 ; i++ ){
     dummy = d[i] - fabs(e[i]) - fabs(e[i+1]);
     shift = min(shift, dummy);
     }
     dummy = d[msize-1] - fabs(e[msize-1]);
     shift = min(shift, dummy);
     
     psgn = 1.0;
     if ( shift < 0.e0){
     psgn = 1.;
     psigma = shift;
     }
     if ( shift == 0.0e0 )
     shift = sqrt(DLAMCHE);
   */
   
   /*
     compute splitting
   */
   
   il = 1;
   iu = 1;
   range = 3;
   order = 1;
   m = 0;
   *info = 0;
   dstebz3_( &range, &order, n, lb, ub, &il, &iu, abstol, d, e+1,
	     &m, nsplit, eval, iblock, isplit, &work[msize], i_work, info);
   
   if ( *info != 0 ) {
     if ( me == 0 )
       printf(" PeIGS error from dstebz %d ...trying dsterf \n", (int)*info );

     file = fopen(filename, "a+");
     fprintf(file, "info = %d \n", (int)linfo);
     fprintf(file, "%d \n", (int)msize);
     for ( iii = 0; iii < msize; iii++)
       fprintf(file, "%d %20.16f %20.16f \n", (int)iii, d[iii], e[iii]);
     fflush(file);
     fclose(file);
     
     for (i = 0; i < msize; i++ )
       work[i] = d[i];
     
     for (i = 0; i < msize; i++ )
       work[msize + i] = e[i];
     
     dsterf_( &msize, &work[0], &work[msize+1], info);
     if ( *info != 0 ){
       if ( me == 0 )
	 printf(" error from dsterf %d \n", (int)*info );
     }
     eval[0] = work[0];
   }
   
   psigma = 0.e0;
   psgn = 1.0;

   /*
   if ( eval[0] < sqrt(sqrt(DLAMCHS)))
     psigma = -fabs(eval[0])-sqrt(sqrt(DLAMCHE));
     */
   
   if ( eval[0] < sqrt(sqrt(DLAMCHU)))
     psigma = -fabs(eval[0])-sqrt(sqrt(DLAMCHE));
   

#ifdef DEBUG99   
   printf("eval %g shift %f sqrt %f  \n", eval[0], psigma, sqrt(sqrt(DLAMCHE)));
#endif
   
   for (i = 0; i < msize; i++ )
     work[i] = d[i] - psgn*psigma;
   
   il = 1;
   iu = msize;
   range = 3;
   order = 1;
   m = 0;
   *info = 0;
   dstebz3_( &range, &order, n, lb, ub, &il, &iu, abstol, work, e+1,
	     &m, nsplit, eval, iblock, isplit, &work[msize], i_work, info);

#ifdef DEBUG99
   for ( jjj = 0; jjj < msize; jjj++ )
     printf(" eval[%d] = %g \n", jjj, eval[jjj]);
   
   for ( jjj = 0; jjj < msize; jjj++ )
     printf(" iblock[%d] = %d \n", jjj, iblock[jjj]);
   for ( jjj = 0; jjj < msize; jjj++ )
     printf(" isplit[%d] = %d \n", jjj, isplit[jjj]);
#endif
   
   if ( *info != 0 ) {
     for (i = 0; i < msize; i++ )
       work[i] = d[i];
     
     for (i = 0; i < msize; i++ )
       work[msize + i] = e[i];
     
     dsterf_( &msize, &work[0], &work[msize+1], info);

     if ( *info != 0 ){
       if ( me == 0 )
	 printf(" error from dsterf %d \n", (int)*info );
     }
     
     for (i = 0; i < msize; i++ ){
       eval[i] = work[i];
     }
   }

#ifdef DEBUG1    
   printf(" psigma = %f psgn = %f onenrm %g  \n", psigma, psgn, onenrm);
#endif
   
   
   /*
     factor shifted matrix into dplus and lplus
   */
   
   /*
    */
   
   for (i = 0; i < *n; i++ )
     work[i] = d[i] - psgn*psigma;
   
   i1split = 0;
   for ( iii = 0; iii < *nsplit; iii++ ){
     jsplit = isplit[iii];
     blksz = jsplit-i1split;
     dptr = &dplus[i1split];
     lptr = &lplus[i1split];
     
     /*
       LDL' factorization of SHIFTED tridiagonal
       */
     
     peigs_tldlfact(&blksz, &work[i1split], &e[i1split], dptr, lptr);
/*
       peigs_dlasq1( blksz, dptr, lptr, &eval[i1split], &work[*n], info );
*/
     
     /*
       for ( jjj = 0; jjj< blksz; jjj++ )
       printf(" jjj = %d dptr %f lptr %f \n", jjj, dptr[jjj], lptr[jjj]);
       */
     
     
     j = iii+1;
/*
     for ( jjj = i1split; jjj < jsplit; jjj++ )
       iblock[jjj] = j;
*/
     i1split = jsplit;
   }
   
   /* The following assumes that the matrix does not split */
   
   return;
}



