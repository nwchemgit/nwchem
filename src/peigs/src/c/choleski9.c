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

void choleski( n, vecA, mapA, iscratch, scratch, info )
     
     Integer            *n, mapA[], iscratch[], *info;
     
     DoublePrecision         scratch[];
     
     DoublePrecision        **vecA;
{
  
  /*
   *  ==================================================================
   *
   *
   *  The original code is due to Mike Heath.  George Fann
   *  bastardized it for our code.  We thank Mike
   *  for his generosity.  Any bugs and problems are introduced
   *  by us.
   * 
   *
   *  Currently, this procedure only handles the case where the
   *  lower triangular part of A is distributed to processors by columns.
   *
   *
   *  Purpose
   *  =======
   *
   *  choleski computes the Cholesky factorization of a real symmetric
   *  positive definite matrix A in packed storage format.
   *
   *  The factorization has the form
   *     A = L  * L'
   *  where L is lower triangular.  On output the matrix A is overwritten by
   *   the Choleski factor L.
   *
   *
   *  Arguments
   *  =========
   *
   *  In the following let:
   *
   *     me     = The id of this processor
   *
   *     nvecsA = The number of columns of A owned by this processor.
   *              In particular, nvecsA = count_list(me, mapA, n).
   *              
   *     nprocsA = Number of distinct processor id's in "mapA".  In
   *               particular, nprocs = reduce_list( n, mapA, iscratch )
   *
   *  n       (input) (Integer *) scaler
   *          The order of the matrix A.  n >= 0.
   *
   *  vecA    (input/output) (DoublePrecision **) array, length (nvecsA)
   *          On entry, the portion of the lower triangular part of the
   *          symmetric matrix A owned by this processor.
   *
   *          vecA[i] points to an array containing the lower triangular
   *          part of the i-th column of A which is owned by this
   *          processor, i = 0 to nvecsA.
   *
   *          On exit, if *info = 0, the factor L from the Cholesky
   *          factorization A = L*L^(T).
   *
   *  mapA    (input) (Integer *) array, length (*n)
   *          mapA[i] = the id of the processor which owns column i of A,
   *          i = 0 to *n - 1.
   *
   *  iscratch  (integer workspace) (Integer *) array, length (3 n + 2)
   *
   *  scratch   (DoublePrecision workspace) (DoublePrecision *) array,
   *            length ( MAX{ n+1, mxlbuf_() / sizeof(DoublePrecision) + 1})
   *
   *  info      (output) (Integer *) scaler
   *            = 0: successful exit
   *            < 0: if info = -k, the k-th argument had an illegal value
   *            > 0: if info = k, 1<= k <= n, the leading minor of order
   *                 k is not positive definite, and the factorization
   *                 could not be completed.
   *
   *
   *
   *  ==================================================================
   */
  
  /*
   *  Local Variables
   */

  char            msg[30];
  Integer             linfo,k, i;
  Integer             me, nprocs, nvecsA;
  Integer            *myvecsA, *iwork;
 
  
  /*
   *  External Procedures
   */
  
  extern Integer      mxmynd_(), mxnprc_();
  extern Integer      count_list(), fil_mapvec_(), reduce_list(), reduce_list2();
  extern void bbcast00();
  
  void     sub_chol0() ;
  extern char    *strcpy();
  
  /*
   *  ---------------------------------------------------------------
   *                      Executable Statements
   *  ---------------------------------------------------------------
   */
  
  /*
   *  Quick Return if possible.
   */
  
  strcpy ( msg, "Error in Choleski\n");
  
  linfo = 0;
  if ( info == NULL ) {
    linfo = -6;
    xerbla_( "Choleski \n", &linfo);
    return;
  }
  
  *info = 0;
  
  if ( n == NULL ) {
    linfo = -1;
    xerbla_( "Choleski \n", &linfo);
    return;
  }
  
 if ( *n < 1) {
    linfo = -1;
    xerbla_( "Choleski \n", &linfo);
    return;
  }
  
  /*
   *  Get this processor's id, and the number of allocated nodes.
   */

  
  me     = mxmynd_();
  nprocs = mxnprc_();
  
  iwork = mapA;
  for ( k = 0; k < *n; k++ ) {
    if ( iwork++ == NULL ) {
      linfo = -2;
      xerbla_ ( "Choleski \n", &linfo);
      return;
    }
  }
  
  /*
   *     Test the input map array
   */
  
  for ( k = 0; k < *n; k++ ) {
    i = mapA[k];
    if (( i  < 0 ) ||   ( i > nprocs - 1 )) {
      linfo = -3;
      xerbla_ ( "Choleski \n", &linfo);
      return;
    }
  }
  
  *iscratch = *n;
  iwork = iscratch + 1;
  for ( k = 0; k < *n; k++ )
    *(iwork++) = mapA[k];
  
  linfo = 0;
  k = (*n + 1)*sizeof(Integer);
  iwork = iscratch + *n + 1;
  pxerbla2_ ( &k, (char *) iscratch, mapA, n, iwork, &linfo );
  
  g_exit_ ( &linfo, "CHOLESKI: Mapping inconsistancies.\n", mapA, n, iwork, scratch);
  
  /*
   *  Partition workspace.
   */

  myvecsA = iscratch;
  nvecsA = fil_mapvec_( &me, n, mapA, myvecsA );
  iwork   = iscratch + nvecsA;
  
  if (nvecsA < 1 ){
    *info = -50;
     return;
  } 
  
  /*
   *  Factor matrix.
   */
  
  sub_chol0( me, *n, vecA, mapA, nvecsA, myvecsA, iwork, scratch, info );
  
  return;
}

  
void sub_chol0( me, n, col, map, ncols, mycols, i_scratch, scratch, info )
     Integer             me, n, ncols;
     Integer            *map, *mycols, *i_scratch, *info;
     DoublePrecision         *scratch;
     DoublePrecision        **col;
{
  
  /*
   *  ==================================================================
   *
   *  Submatrix-Cholesky factorization.
   *  Uses fan-out communication without send-ahead.
   *
   *  The original code is due to Mike Heath.  George Fann
   *  bastardized it for our code.
   *
   *  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   *  This procedure does no input error checking and should only be
   *  called via "choleski".
   *  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   *
   *  Workspace requirements:
   *  -----------------------
   *
   *     Variable       Type       Length
   *    -----------  ----------    ------
   *    i_scratch     (Integer *)   2*n+1
   *
   *    scratch       (DoublePrecision *)    n+1
   *
   *
   *  On exit:
   *  --------
   *
   *    *info = 0  means factorization finished successfully
   *    *info = k  1 <= k <= n, means the leading minor of order k is
   *               not positive definite, and the factorization could not
   *               be completed.
   *
   *  ==================================================================
   */
  
  /*
    *  Local Variables
    */

   static Integer      IONE = 1;
   static DoublePrecision   ZERO = 0.0e0, ONE = 1.0e0;

   Integer             i, j, k, m, indx, num_procs;
   DoublePrecision          t, *q;
   long length, root;

   /*
    *  External Procedures
    */

   extern void     dscal_(), daxpy_(), dcopy_();
   extern void     chol_pipe_bcast();
   extern Integer  reduce_list2();
   extern void bbcast00();

/*
   extern void iCC_bcast(), iCC_work_init();
*/
   
   /*
    *  Intrinsic Procedures
    */

   /*
    *  ---------------------------------------------------------------
    *                      Executable Statements
    *  ---------------------------------------------------------------
    */

  *info = 0;


#ifdef DEBUG1
   for ( j = 0; j < n; j++ )
     fprintf(stderr, " me = %d j = %d map %d \n", me, j, map[j]);
#endif


   num_procs = reduce_list2( n, map, i_scratch );

#ifdef DEBUG1
   for ( j = 0; j < num_procs; j++ )
     fprintf(stderr, " me = %d j = %d map %d num_procs = %d \n", me, j,
	     map[j], num_procs);
#endif
   

   j = 0;   

   /*
     #ifdef Intel
     i_random(num_procs, i_scratch, i_scratch+num_procs+1);
     #endif
     */
   
   /*
     iCC_work_init(( unsigned long ) n);
     */
   
   for ( k = 0; k < n; k++ )
     {
       m = n - k;
       if ( map[ k ] == me )
	 {
	   /*
	    *  Check for postive definiteness and scale column k.
	    */
	   
	   if ( *col[ j ] <= ZERO ){
	     scratch[ m ] = ONE;
	   }
	   else
	     {
	       t = ONE/sqrt( *col[ j ] );
	       dscal_( &m, &t, col[ j ], &IONE );
	       dcopy_( &m, col[ j ], &IONE, scratch, &IONE );
	       j++;
	       scratch[ m ] = ZERO;
	     }
	 }
       
       /*
	*  Distribute new column k, and info.
	*/
       
       i = (m+1) * sizeof(DoublePrecision);
       length = (long) i;
       root = (long) map[k];
       
       if ( k > -1 ) {
	 /*
	   iCC_bcast( scratch, length, root );
	   */
	 bbcast00((char *) scratch, i, k, map[k], num_procs, i_scratch );
       }
       else
	 chol_pipe_bcast( (char *) scratch, i, k, map[ k ], n - k, &map[ k ], i_scratch );
       
       if ( scratch[ m ] != ZERO )
	 {
	   *info = k + 1;
	   break;
	 }
       
       /*
	*  Update the rest of the matrix.
	*/
       
       for ( i = j; i < ncols; i++ )
	 {
	   indx = mycols[ i ];
	   q    = scratch + indx - k;
	   m    = n - indx;
	   t    = -(*q);
	   daxpy_( &m, &t, q, &IONE, col[ i ], &IONE );
	 }
     }
   
   /*
    *  Global check of info.  Just pipeline info from owner of column,
    *  n-1, who always knows the correct final value of info, to the
    *  rest of the processors in "map".
    */
   
   /*
     chol_pipe_bcast( info, sizeof(Integer), n+1, map[ n-1 ], n, map, i_scratch );
     */
   
   num_procs = reduce_list2( n, map, i_scratch );   
#ifdef MESH
   i_random(num_procs, i_scratch, i_scratch+num_procs+1);   
#endif
   bbcast00( (char *) info, sizeof(Integer), n+1, map[n-1], num_procs, i_scratch );
   
   /*
     length = (long) sizeof(Integer);
     root = (long) map[0];
     iCC_bcast( info, length, root );
     */
   
   return;
 }


