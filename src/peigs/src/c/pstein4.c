
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
 *  -- PEIGS  routine (version 3.1) --
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

/*
  static double peigs_timer;
  */


void pstein4 ( n, dd, ee, dplus, lplus, ld, lld, meigval, eval, iblock, nsplit, isplit,
	       mapZ, vecZ, clustr_info, ddwork, iiwork, ppiwork, info )
     
     Integer            *n, *meigval, *iblock, *nsplit, *isplit,
  *mapZ, *clustr_info,
  *iiwork, *info, **ppiwork;
     DoublePrecision        *dd, *ee, *dplus, *lplus, *ld, *lld, *eval, **vecZ, *ddwork;
{
  
  /*
   *
   * Driver to compute eigenvectors of a symmetric tridiagonal matrix.
   * Basically a parallel version of LAPACK's DSTEIN.
   *
   *  n ......... The dimension of the matrix.
   *
   *  dd ........ The diagonal of the matrix.
   *
   *  ee ........ The off-diagonal of the matrix.  ee[0] is junk, the
   *              actual off-diagonal starts at ee[1].
   *
   *  meigval ... (input) Integer
   *              The number of eigenvalues found
   *
   *  eval ...... (input) DoublePrecision array, length (neigval)
   *              The actual eigenvalues.  Sorted by split of block,
   *              and within a block sorted by value.
   *
   *  iblock .... (input) Integer array, length (neigval)
   *              Array from pstebz.
   *
   *  nsplit .... (input) Integer
   *              The number of split points in the tridiagonal matrix.
   *              From pstebz.
   *
   *  isplit .... (input) Integer array, length (nsplit)
   *              Array from pstebz.
   *
   *  mapZ    (input/output) INTEGER array, dimension (meigval)
   *          On entry:
   *
   *          mapZ(i) = the id of a processor which has room for the
   *                    i-th eigenvector, i = 0 to meigval-1
   *          On exit:
   *            mapZ(i) = the id of the processor which actually owns the i-th
   *                      eigenvector, i = 0 to meigval-1.
   *
   *          The value of mapZ on exit may be different then on entry.
   *
   *  vecZ    (output) array of pointers to DoublePrecision (DoublePrecision **)
   *                   dimension ( nvecsZ )
   *          On entry:
   *             vecZ[i], i = 0 to nvecsZ-1, should point to an array of length n.
   *
   *          On exit:
   *
   *            vecZ[i], i = 0 to nvecsZ-1, points to the i-th eigenvector
   *            (as determined by the exit values in mapZ) owned by this
   *            processor.
   *
   *            The eigenvectors are normalized such that: Z'* Z = I.
   *
   *
   *  ddwork .... (workspace) DoublePrecision array
   *  iiwork .... (workspace) Integer array,
   *  ppiwork ... (workspace) array of pointers to Integer,
   *
   *
   *  INFO ... (output) INTEGER
   *
   *           PSTEIN attempts to return the same INFO on all processors in MAPZ.
   *           Currently, however, if the input data is invalid, -50 < INFO < 0,
   *           then INFO will be different on different processors.
   *
   *           = 0:   successful exit
   *
   *           -50 <= INFO < 0:
   *                  Then the (-INFO)-th argument had an illegal value
   *
   *           = -51: Then the input data was not the same on all processors in MAPZ
   *
   *           0 < INFO <= MEIGVAL:
   *                  Then the INFO-th eigenvector failed to converge, but all
   *                  lower numbered eigenvectors did converge.
   *
   *           N < INFO:
   *                  Then the residual,
   *                  res == ( max_i |T z_i-lamba_i z_i | /( eps |T|) ),
   *                  for the tridiagonal eigenproblem was excessively large.
   *                  In particular, INFO = N + (Integer) log10( res )
   */
  
  /*
  static Integer      IONE = 1;
  static DoublePrecision   ZERO = 0.0e0;
  */
  
  Integer              k, ii, isize, me, msize,
                   nproc, nn_proc, neigval, nvecsZ,
                   ibad, linfo, maxinfo, imin, nacluster;
  
   Integer            *iwork, *proclist,
                  *i_scrat, *iscrat, *mapvZ,
                  *iscratch, *icsplit;

   char            msg[ 25 ];
   char            msg2[ 25 ];

   Integer           **piwork, max_sz, sync_proc, numclstr;

   DoublePrecision         *dwork, *d_scrat, dbad[1], syncco[1];
   extern DoublePrecision tcgtime_();


  extern Integer      mxnprc_(), mxmynd_();

  extern Integer clustrf4_();
  
  extern Integer      reduce_list2(), count_list(), clustrinv4_();
  extern Integer      indxL ();

  extern void reduce_maps();
  extern void     xstop_();
  extern void     pgexit();
  extern void     pdiff();
  extern void gsum01(), gmax00();
  extern void bbcast00(), tresid();
  extern Integer       fil_mapvec_(), clustrfix_();
  
/*
 *  ---------------------------------------------------------------
 *                      Executable Statements
 *  ---------------------------------------------------------------
 */

   /*
    *  Get this processor's id, and the number of allocated nodes.
    */

  me = mxmynd_();
  nproc = mxnprc_();

  /*
  peigs_timer = tcgtime_();
  */
  
  strcpy( msg,  "Error in pstein4." );
  
  
#ifdef DEBUG1
   fprintf(stderr, "me = %d In pstein \n", me );
#endif

   /*
    *     Test the input parameters.
    */
  
  linfo = 0;
  
  if ( n == NULL )
    linfo = -1;
  else if ( dd == NULL )
    linfo = -2;
  else if ( ee == NULL )
    linfo = -3;
  else if ( meigval == NULL )
    linfo = -4;
  else if ( eval == NULL )
    linfo = -5;
  else if ( iblock == NULL )
    linfo = -6;
  else if ( nsplit == NULL )
    linfo = -7;
  else if ( isplit == NULL )
    linfo = -8;
  else if ( mapZ == NULL )
    linfo = -9;
  else if ( vecZ == NULL )
    linfo = -10;
  else if ( ddwork == NULL )
    linfo = -11;
  else if ( iiwork == NULL )
    linfo = -12;
  else if ( ppiwork == NULL )
    linfo = -13;
  else if ( info == NULL )
    linfo = -14;
  
  if ( linfo != 0 ) {
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
  
  if ( *n == 0  ||  *meigval == 0 )
    return;
  
  /*
   *  Continue error checking.
   */
  
  if ( *n < 0 )
    *info = -1;
  else if ( *meigval < 0  ||  *meigval > *n )
    *info = -4;
  else if ( iblock[ 0 ] < 1  ||  iblock[ *meigval - 1 ] > *nsplit )
    *info = -6;
  else if ( *nsplit < 1  || *nsplit > *n )
    *info = -7;
  else if ( isplit[ 0 ] < 1  ||  isplit[ *nsplit - 1 ] != *n )
    *info = -8;
  
  if ( *info == 0 )
    for ( k = 1; k < *meigval; k++ )
      if ( iblock[ k ] == iblock[ k - 1 ]  &&  eval[ k ] < eval[ k - 1 ] )
	*info = -5;
  
  if ( *info == 0 )
    for ( k = 1; k < *meigval; k++ )
      if ( iblock[ k ] < iblock[ k - 1 ] )
	*info = -6;
  
  if ( *info == 0 )
    for ( k = 1; k < *nsplit; k++ )
      if ( isplit[ k ] <= isplit[ k - 1 ] )
	*info = -8;
  
  if ( *info == 0 )
    for ( k = 0; k < *meigval; k++ )
      if ( mapZ[ k ] < 0  ||  mapZ[ k ] > nproc - 1 )
	*info = -9;
  
  if ( *info == 0 ) {
      /*
       * Count the number of columns of Z owned by this processor.
       * Must own something.
       */
      
      nvecsZ = count_list( me, mapZ, meigval );
      if ( nvecsZ <= 0 )
	return;

      for ( k = 0; k < nvecsZ; k++ )
        if ( vecZ[ k ] == NULL )
	  *info = -10;
  }
  
  if( *info != 0 ) {
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
  
  proclist = iiwork;
  reduce_maps( *meigval, mapZ, 0, mapZ, 0, mapZ, &nn_proc, proclist );
  iscratch = iiwork + nn_proc;
  iscratch += 4 * msize;

  /*
   *  Check scaler inputs.
   */
  
  iscratch [ 0 ] = *n;
  iscratch [ 1 ] = *meigval;
  iscratch [ 2 ] = *nsplit;
  
  isize = 3 * sizeof( Integer );
  strcpy( msg2, "n,meigval,or nsplit " );
  pdiff( &isize, (char *) iscratch, proclist, &nn_proc, iscratch+3, msg, msg2, &linfo );
  
  pgexit( &linfo, msg, proclist, &nn_proc, ddwork);
  
  if ( linfo != 0 ) {
    *info = -51;
    return;
  }
  
  /*
   *  Check rest of inputs.
   */
  
  maxinfo = 0;
  
  isize = msize * sizeof( DoublePrecision );
  strcpy(msg2, "dd ");
  pdiff( &isize, (char *) dd, proclist, &nn_proc, (Integer *) ddwork, msg, msg2 , &linfo );
  maxinfo = max( maxinfo, linfo );
  
  strcpy(msg2, "ee ");
  pdiff( &isize, (char *) ee, proclist, &nn_proc, (Integer *) ddwork, msg, msg2, &linfo );
  maxinfo = max( maxinfo, linfo );
  
  strcpy(msg2, "eval ");
  isize = *meigval * sizeof( DoublePrecision );
  pdiff( &isize, (char *) eval, proclist, &nn_proc, (Integer *) ddwork, msg, msg2, &linfo );
  maxinfo = max( maxinfo, linfo );
  
  isize = *meigval * sizeof( Integer );
  strcpy(msg2, "iblock ");  
  pdiff( &isize, (char *) iblock, proclist, &nn_proc, iscratch, msg, msg2, &linfo );
  maxinfo = max( maxinfo, linfo );
  
  isize = *nsplit * sizeof( Integer );
  strcpy(msg2, "isplit ");  
  pdiff( &isize, (char *) isplit, proclist, &nn_proc, iscratch, msg, msg2, &linfo );
  maxinfo = max( maxinfo, linfo );
  
  isize = *meigval * sizeof( Integer );
  strcpy(msg2, "mapZ ");  
  pdiff( &isize, (char *) mapZ, proclist, &nn_proc, iscratch, msg, msg2, &linfo );
  maxinfo = max( maxinfo, linfo );
  
  linfo = maxinfo;
  
  pgexit( &linfo, msg, proclist, &nn_proc, ddwork);
  
  if ( linfo != 0 ) {
    *info = -51;
      return;
  }
  
  /*
    * ----------------------------------------
    * All input data is good. Start computing. 
    * ----------------------------------------
    */
  
  neigval = *meigval;

/*
  iwork  = iiwork;
*/
  iwork  = iscratch;
  dwork  = ddwork;
  piwork = ppiwork;
  
  /*
   * Set up proclist work array.
   */
  
  nn_proc = reduce_list2( neigval, mapZ, iwork );
  proclist = iwork;
  iscrat =   iwork + nn_proc ;
  
  /*
   * Set up remaining integer work arrays.
   */
  
  iwork += (6 * msize + 1);
  
  /* iscrat is not used, so take it over for icsplit which is of length neigval. */
  
  icsplit = iscrat;
  iscrat += 6*msize + 1;
  
  
  mapvZ   = iwork;
  iwork  += nvecsZ;
  i_scrat = iwork;
  
  /*
   * Set up DoublePrecision precision work arrays.
   */
  
  d_scrat = dwork;
  
  /*
   * dscrat is of length 5 * msize
   */
  
  /*
   * Set up pointer to integer work arrays.
   */
  
  nacluster = 0;
  numclstr = 0;
  isize = 0;
  
  /*
    isize = clustrf4_(&msize, dplus, lplus, &neigval,
    eval, mapZ, vecZ, iblock, nsplit, isplit,
    clustr_info,
    &numclstr, icsplit, i_scrat);
    */
  
  isize = clustrfix_(msize, dd, ee, neigval, eval, iblock, *nsplit, isplit,
		     &numclstr, clustr_info );
  
  
  /*
    for ( ii = 0; ii < 4*numclstr; ii++ ){
    printf(" pstein4 me = %d ii %d clustrf_info %d \n", me, ii, clustr_info[ii]);
    }
    */
  
  
  nacluster = numclstr;
  
  
  if ( isize < 0 ){
    *info = -99;
    fprintf(stderr, " Node %d: error in clustrf isize = %d neigval = %d \n", 
	    me, isize, neigval );
    xstop_( info );
  }
  
  nvecsZ = fil_mapvec_( &me, &neigval, mapZ, mapvZ );
  
  /*
   * syncronize processors
   */
  
  max_sz = 0;
  for ( ii = 0; ii < nproc; ii++ )
    max_sz = max( i_scrat[ii], max_sz );
  
  if( max_sz == 0 ) {
    sync_proc = proclist[0];
  }
  else {
    sync_proc = -1;
    for ( ii = 0; ii < nproc; ii++ ) {
      if ( i_scrat[ii] == max_sz ) 
        sync_proc = max(sync_proc, ii);
    }
  }
  
  /*
    bbcast00( (char *) d_scrat, 1, 999, sync_proc, nn_proc, proclist);
  */
  
  
  syncco[0] = 0.0e0;
  gsum00( (char *) syncco, 1, 5, 10, mapZ[0], nn_proc, proclist, d_scrat);

  
  /*
   * Compute eigenvectors
   */
  
#ifdef DEBUG  
  printf(" in pstein4 me = %d \n", me);
#endif

  /*
    printf(" in pstein4 bbcast00 me = %d numclstr %d \n", me, numclstr);
    fflush(stdout);
  */
  
  nacluster= numclstr;
  
  /*
    bbcast00( (char *) &clustr_info[0], 4*msize*sizeof(Integer), 9, proclist[0], nn_proc, proclist);
  */
  
  /*
    printf(" in pstein4 after bbcast00 me = %d \n", me);
    fflush(stdout);
  */
  
  ibad = 0;
  ibad = clustrinv4_( &msize, dd, ee, dplus, lplus, ld, lld,
		      eval, clustr_info, &numclstr, mapZ,
		      mapvZ, vecZ, &imin, &nacluster, icsplit,
		      i_scrat, d_scrat);
  
  /*
   * syncronize processors
   */
  
  /*
   * Get same value of ibad for all processors in proclist.
   */
  
	  

  
  dbad[0] = (DoublePrecision) ibad;
  gmax00( (char *) dbad, 1, 1, 888, proclist[0], nn_proc, proclist, dwork );

  ibad = (Integer) dbad[0];
  ibad = -ibad;
    
  if( ibad == msize+10 )
    ibad = 0;
  
  /*
   * Check residual of tridiagonal eigenproblem.
   */
  
  
  *info = ibad;
  
  
#ifdef DEBUG7
  printf(" out pstein4 me = %d \n", me);
#endif
  
  return;
}






