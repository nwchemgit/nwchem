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
#include <stdlib.h>

#include "globalp.c.h"

#define MSG_START 25000
#define ZERO ((DoublePrecision) 0.0e0)

#define max(a,b) ((a) > (b) ? (a) : (b))

void residual2( n, colA, mapA, colB, mapB, m, colZ, mapZ, eval,
                ibuffptr, iwork, work, res, info)
     Integer *n, *mapA, *mapB, *m, *mapZ, *iwork, *info;
     DoublePrecision **colA, **colB, **colZ, *eval, **ibuffptr, *work, *res;
     
/*
       
   This subroutine computes the residual
       
   res = max_{i} | A z_{i} - \lambda_{i} B z_{i} |/( | A | * ulp )
       
   where 

   A is an n-by-n symmetric matrix, in packed storage format,
   column (or equivalently row) distributed
   
   B is an n-by-n symmetric matrix, in packed storage format,
   column (or equivalently row) distributed
   
   (lambda_i, z_i) is a generalized eigen-pair and
   Z is an n-by-m matrix of eigenvectors, in packed storage format,
   column distributed
   
   ULP = (machine precision) * (machine base)
       
   |A z_{i} ... | is the infinity-norm,
   |A| is the 1-norm of A,

   res is reasonable if it is of order about 50 or less.
       
   Arguments
   ---------
   In the following:

     INTEGER          = "pointer to Integer"
     DOUBLE PRECISION = "pointer to DoublePrecision"

     me     = this processor's id (= mxmynd_())
     nprocs = number of allocated processors ( = mxnprc_())
     nvecsA = number of entries in mapA equal to me
              (= count_list( me, mapA, n ))
     nvecsB = number of entries in mapB equal to me
              (= count_list( me, mapB, n ))
     nvecsZ = number of entries in mapZ equal to me
              (= count_list( me, mapZ, m ))
     sDP    = sizeof( DoublePrecision )

       
   n....... (input) INTEGER
            size of the symmetric matrices A and B,
            and the number of rows in Z.

   colA ... (input) array of pointers to DoublePrecision,
                    length (nvecsA)
            The part of matrix A owned by this processer stored
            in packed format, i.e., colA[i] points to the diagonal
            element of the i-th column (or equivalently row) of A
            owned by this processor, i = 0 to nvecsA-1.
                
   mapA ... (input) INTEGER array, length (n)
            The i-th column (or equivalently row) of A is 
            owned by processor mapA[i], i = 0 to n-1.

   colB ... (input) array of pointers to DoublePrecision,
                    length (nvecsB)
            The part of matrix B owned by this processer stored
            in packed format, i.e., colB[i] points to the diagonal
            element of the i-th column (or equivalently row) of B
            owned by this processor, i = 0 to nvecsB-1.
                
   mapB ... (input) INTEGER array, length (n)
            The i-th column (or equivalently row) of B is 
            owned by processor mapB[i], i = 0 to n-1.

   m....... (input) INTEGER
            number of columns of Z, i.e., # of eigenvalues/eigenvectors

   colZ ... (input) array of pointers to DoublePrecision,
                    length (nvecsZ)
            The part of matrix Z owned by this processer, stored
            in packed format, i.e., colZ[i] points to the start
            of the i-th column of matrix Z owned by this
            processor, i = 0 to nvecsZ-1.
                
   mapZ ... (input) INTEGER array, length (m)
            The i-th column of matrix Z is owned by processor
            mapZ[i], i = 0 to m-1.

   eval.... (input) DOUBLE PRECISION array, length (m)
            the eigenvalues
       
   ibuffptr (workspace) array of pointers to DoublePrecision,
                        length (3 * nvecsZ + 3)

   iwork... (workspace) INTEGER array, 
                        length( m + nvecsA+nvecsB+nvecsZ+
                                MAX{ i_mxm88, nprocs, nvecsA + n )
                        where i_mxm88 = 2*(n+m)+nvecsA+nvecsZ+nvecsB

   work.... (workspace) DOUBLE PRECISION array,
                        length( 2*nvecsZ * n + maximum( d_mxm88,
                                                      n + 1 + mxlbuf_() / sDP + 1 )
                        where d_mxm88 = maximum ( (nvecsZ+1)*n, mxlbuf_()/sDP + 1 )
       
   res..... (output) INTEGER
            the residual described above.

   info.... (output) INTEGER
            = 0, not currently used
       
 */
{
  
  static Integer IONE = 1;
  
  Integer ll, i, *iscrat, *mapvecA, *mapvecB, *mapvecZ;
  Integer nvecsA, nvecsB, nvecsZ;
  Integer me;
  Integer nprocs, *proclist;
  Integer k;
  
  DoublePrecision t, derror;
  DoublePrecision normA, ulp;
  
  DoublePrecision *ptr, *scrat;
  
  /*
    copies of the vecZ matrix
    */
  
  DoublePrecision **dbuffptr, **vecZ1, **vecZ2;
  
  extern void gmax00();
  
  
  /*
    blas call
    */
  
  extern DoublePrecision dnrm2_();
  extern DoublePrecision ddot_();
  extern void dcopy_();
  extern void daxpy_();
  extern DoublePrecision damax_();
  
  extern Integer mxwrit_(), mxread_();
  extern Integer menode_(), mxmynd_();
  
  extern void chol_pipe_bcast();
  extern Integer fil_mapvec_();
  extern Integer count_list();
  extern void mxm88();
  extern void sonenrm ();
  extern Integer reduce_list2();
  extern void    mdiff1_(), mdiff2_(), bbcast00();


  /*
    usual story about error handling
    */
  
  ll = *n;
  
  if ( ll < 1 )
    return;  /* error */
  
  /*
    should perform global sync to check for errors
    */
  
  me = mxmynd_();
  
  iscrat = iwork;  
  mapvecA = iscrat;
  nvecsA = fil_mapvec_( &me, &ll, mapA, mapvecA);
  iscrat += nvecsA;
  
  mapvecB = iscrat;
  nvecsB = fil_mapvec_( &me, &ll, mapB, mapvecB);
  iscrat += nvecsB;
  
  mapvecZ = iscrat;
  nvecsZ = fil_mapvec_( &me, m, mapZ, mapvecZ );
  iscrat += nvecsZ;
  
  proclist = iscrat;
  nprocs = reduce_list2( *m, mapZ, proclist );
  iscrat += nprocs;
  
  if (( nvecsA == 0 ) && ( nvecsB == 0 ) && ( nvecsZ == 0 ))
    return;
  
  scrat = work;
  
  vecZ1 = (DoublePrecision **) malloc( nvecsZ*sizeof(DoublePrecision *));
  vecZ2 = (DoublePrecision **) malloc( nvecsZ*sizeof(DoublePrecision *));
  
  dbuffptr = ibuffptr;
  
  /*
    vecZ1 = dbuffptr;
    dbuffptr += nvecsZ + 1;

    vecZ2 = dbuffptr;
    dbuffptr += nvecsZ + 1;
    */
  
  /*
    copy over the matrix to a buffer area
    */
  
  ptr = scrat;
  for ( i = 0; i < nvecsZ; i++ ) {
    vecZ1[i] = ptr;
    dcopy_( &ll, colZ[i], &IONE, ptr, &IONE );
    ptr += ll;
  }
  
  for ( i = 0; i < nvecsZ; i++ ) {
    vecZ2[i] = ptr;
    dcopy_( &ll, colZ[i], &IONE, ptr, &IONE );
    ptr += ll;
  }
  
  scrat = ptr;
  
  /*
    compute vecZ1 <- A . Z
    */
  
  mxm88( &ll, colA, mapA, m, vecZ1, mapZ, iscrat, scrat, dbuffptr);
  
  /*
    compute vecZ2 <- B . Z
    */
  
  mxm88( &ll, colB, mapB, m, vecZ2, mapZ, iscrat, scrat, dbuffptr);
  
  k = 0;
  for ( i = 0; i < *m; i++ ) {
    if ( mapZ[i] == me ) {
      t = -eval[i];
      daxpy_( &ll, &t, vecZ2[k], &IONE, vecZ1[k], &IONE);
      k++;
    }
  }
  
  /*
    vecZ1[i] is now A z_{i} - lambda_{i} B z_{i}
    */
  
  derror =  0.0e0;
  for ( i = 0; i < nvecsZ; i++ ) {
    derror = max( damax_( &ll, vecZ1[i], &IONE), derror);
  }

  gmax00( (char *) &derror, 1, 5, 16, proclist[0], nprocs, proclist, scrat);
  
  sonenrm( &ll, colA, mapA, &normA, iscrat, scrat, info);
  
  if( normA == 0.0e0 ) normA = 1.0e0;

  /*
   * Make sure all processors in mapZ know normA.
   */

  mdiff1_( m, mapZ, n, mapA, iscrat, &nprocs ); 

  if( nprocs > 0 ){
      iscrat[nprocs] = mapA[0];
      nprocs++;
      bbcast00( (char *) &normA, sizeof(DoublePrecision), 2, mapA[0], nprocs, iscrat );
  }
    
  if( nvecsZ > 0 ) {

    ulp = DLAMCHE * DLAMCHB ;

    *res = derror / normA / ulp;

  }
  
  /*
   * Make sure all processors in mapA, mapB and mapZ know res.
   */

  mdiff2_( n, mapA, n, mapB, m, mapZ, iscrat, &nprocs ); 

  if( nprocs > 0 ){
      iscrat[nprocs] = mapZ[0];
      nprocs++;
      bbcast00( (char *) res, sizeof(DoublePrecision), 3, mapZ[0], nprocs, iscrat );
  }

  free(vecZ2);
  free(vecZ1);

  return;
}

