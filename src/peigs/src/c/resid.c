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

#include "globalp.c.h"

#define ZERO ((DoublePrecision) 0.0e0)

#define max(a,b) ((a) > (b) ? (a) : (b))

void resid( n, colA, mapA, m, colZ, mapZ, eval, ibuffptr, iwork, work, res, info)
     Integer *n, *mapA, *m, *mapZ, *iwork, *info;
     DoublePrecision **colA, **colZ, *eval, **ibuffptr, *work, *res;
     
/*
   this subroutine computes the residual
       
   res = max_{i} | A z_{i} - \lambda_{i} z_{i} |/( | A | * ulp )
       
   where 
   A is an n-by-n symmetric matrix, in packed storage format,
   column (or equivalently row) distributed

   (lambda_i, z_i) is a standard eigen-pair of A and
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
     nvecsZ = number of entries in mapZ equal to me
                  (= count_list( me, mapZ, m ))
     sDP    = sizeof( DoublePrecision )

       
   n....... (input) INTEGER
            dimension of the matrix A

   colA ... (input) array of pointers to DoublePrecision,
                    length (nvecsA)
            The part of matrix A owned by this processer stored
            in packed format, i.e., colA[i] points to the diagonal
            element of the i-th column (or equivalently row) of A
            owned by this processor, i = 0 to nvecsA-1.
                
   mapA ... (input) INTEGER array, length (n)
            The i-th column (or equivalently row) of A is 
            owned by processor mapA[i], i = 0 to n-1.

   m....... (input) INTEGER
            number of columns of Z, i.e. # of eigenvalues/eigenvectors

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
                        length (2 * nvecsZ + 1)

   iwork... (workspace) INTEGER array, 
                        length( m + maximum( nprocs, n+nvecsA, i_mxm88 ))
                        where i_mxm88 = 2*(n+m)+nvecsA+nvecsZ

   work.... (workspace) DOUBLE PRECISION array,
                        length( nvecsZ * n + maximum( d_mxm88,
                                                      n + 1 + mxlbuf_() / sDP + 1 )
                        where d_mxm88 = maximum ( (nvecsZ+1)*n, mxlbuf_()/sDP + 1 )
       
   res..... (output) DOUBLE PRECISION
            the residual described above.

   info.... (output) INTEGER
            = 0, not currently used
       
 */
{
  static Integer IONE = 1;
  
  Integer           ll, i, nvecsA, nvecsZ, me, nprocs, k;
  Integer          *iscrat, *proclist;
  DoublePrecision   t, derror, normA, ulp;
  DoublePrecision  *ptr, *scrat;
  DoublePrecision **dbuffptr, **vecZ1;
  
  extern void             dcopy_(), daxpy_();
  extern DoublePrecision  damax_();
  extern Integer          mxmynd_();
  extern Integer          count_list(), reduce_list2();
  extern void             gmax00();
  extern void             mxm88(), sonenrm ();
  extern void             mdiff1_(), bbcast00();

  /*
    usual story about error handling
    */
  
  ll = *n;
  
  *res = 0.e0;
  *info = 0;

  if ( ll < 1 )
    return;  /* error */
  
  /*
    should perform global sync to check for errors
    */
  
  me = mxmynd_();
  
  iscrat = iwork;  

  nvecsA = count_list( me, mapA, &ll );
  nvecsZ = count_list( me, mapZ, m   );
  
  if (( nvecsA == 0 ) && ( nvecsZ == 0 ))
    return;
  
  proclist = iscrat;
  nprocs = reduce_list2( *m, mapZ, proclist );
  iscrat += nprocs;
  
  scrat = work;

  dbuffptr = ibuffptr;

  vecZ1 = dbuffptr;
  dbuffptr += nvecsZ + 1;

  /*
    copy over the matrix to a buffer area
    */
  
  ptr = scrat;
  for ( i = 0; i < nvecsZ; i++ ) {
    vecZ1[i] = ptr;
    dcopy_( &ll, colZ[i], &IONE, ptr, &IONE );
    ptr += ll;
  }
  
  scrat = ptr;
  
  /*
    compute vecZ1 <- A . Z
    */
  
  mxm88( &ll, colA, mapA, m, vecZ1, mapZ, iscrat, scrat, dbuffptr);
  
  k = 0;
  for ( i = 0; i < *m; i++ ) {
    if ( mapZ[i] == me ) {
      t = -eval[i];
      daxpy_( &ll, &t, colZ[k], &IONE, vecZ1[k], &IONE);
      k++;
    }
  }
  
  /*
   * vecZ1[i] is now A z_{i} - lambda_{i} z_{i}
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
  
  if( nvecsZ > 0 ){
    ulp = DLAMCHE * DLAMCHB ;
    *res = derror / normA / ulp;
  }
  
  /*
   * Make sure all processors in mapA and mapZ know res.
   */
  
  mdiff1_( n, mapA, m, mapZ, iscrat, &nprocs ); 
  
  if( nprocs > 0 ){
    iscrat[nprocs] = mapZ[0];
    nprocs++;
    bbcast00( (char *) res, sizeof(DoublePrecision), 3, mapZ[0], nprocs, iscrat );
  }
  return;
}

