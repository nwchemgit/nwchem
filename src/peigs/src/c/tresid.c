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

#define ZERO ((DoublePrecision) 0.0e0)

#define max(a,b) ((a) > (b) ? (a) : (b))

void tresid( n, m, d, e, colZ, mapZ, eval, iwork, work, res, info)
     Integer *n, *m, *mapZ, *iwork, *info;
     DoublePrecision *d, *e, **colZ, *eval, *work, *res;
     
/*
       this subroutine computes the residual
       
       res = max_{i} | T z_{i} - \lambda_{i} z_{i} |/( | T | * ulp )
       
       where T is an n-by-n  tridiagonal matrix,
        ( \lambda_{i} , z_{i} ) is a standard eigen-pair, and
       ULP = (machine precision) * (machine base)
       
       |T| is the 1-norm of T
       |T z_{i} .... | is the infinity-norm
       res is reasonable if it is of order about 50 or less.
       
   Arguments
   ---------
   In the following:

     INTEGER          = "pointer to Integer"
     DOUBLE PRECISION = "pointer to DoublePrecision"

     me     = this processor's id (= mxmynd_())
     nprocs = number of allocated processors ( = mxnprc_())
     nvecsZ = number of entries in mapZ equal to me
                  (= count_list( me, mapZ, n ))
     sDP    = sizeof( DoublePrecision )

       
   n....... (input) INTEGER
            dimension of the matrix T

   m....... (input) INTEGER
            number of eigenvalues/eigenvectors

   d....... (input) DOUBLE PRECISION array, length (n)
            diagonal of T

   e....... (input) DOUBLE PRECISION array, length (n)
            e[0] = junk,
            e[1:n-1] = sub-diagonal of T
                     = super-diagonal of T
  
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
       
   iwork... (workspace) INTEGER array, length(m)

   work.... (workspace) DOUBLE PRECISION array,
                        length( mxlbuf_() / sDP + 1 )
       
   res..... (output) INTEGER
            the residual described above.

   info.... (output) INTEGER
            = 0, not currently used
*/

{
  
  Integer           ll, i, j, nvecsZ, me, nprocs, k;
  Integer          *iscrat, *proclist;
  DoublePrecision   t, derror, normA, ulp, dumb;
  DoublePrecision  *ptr, *scrat;
  
  extern Integer          mxmynd_();
  extern Integer          count_list(), reduce_list2();
  extern void             gmax00();

  
  /*
    usual story about error handling
    should perform global sync to check for errors
    */
  
  ll = *n;

  *info = 0;
  *res = 0.e0;
  
  if ( ll < 1 )
    return;

  me = mxmynd_();
  
  iscrat = iwork;  
  scrat  = work;

  nvecsZ = count_list( me, mapZ, m   );
  
  if( nvecsZ == 0 )
    return;
  
  proclist = iscrat;
  nprocs = reduce_list2( *m, mapZ, proclist );

  if( ll > 1 ) {
    normA = fabs( d[0] ) + fabs( e[1] );
    normA = max(normA, fabs(d[*n-1]) + fabs(e[*n-1]) );
    for (i = 1; i < *n-1; ++i)
      normA = max( normA, fabs(d[i]) + fabs(e[i]) + fabs(e[i+1]) );
    
  }
  else {
    normA = fabs( d[0] );
  }
  if( normA == 0.0e0 ) normA = 1.0e0;

  ulp = DLAMCHE * DLAMCHB ;  


  k      = 0;
  derror =  0.0e0;
  
  for ( i = 0; i < *m; i++ ) {
    if ( mapZ[i] == me ) {
      ptr = colZ[k];
      if( ll > 1 ) {
        t = fabs( ( d[0] * ptr[0] + e[1] * ptr[1] ) - eval[i] * ptr[0] );
        t = max( t, fabs( ( d[ll-1] * ptr[ll-1] + e[ll-1] * ptr[ll-2] ) - eval[i] * ptr[ll-1] ));
        for ( j = 1; j < ll-1; j++ )
	  t = max( t, fabs( ( e[j] * ptr[j-1] + d[j] * ptr[j] + e[j+1] * ptr[j+1] ) - eval[i] * ptr[j] ));
      } 
      else {
        t  = fabs( d[0] * ptr[0] - eval[i] * ptr[0] );
      } 

      dumb = t/normA/ulp;
      derror = max( t, derror);
      k++;
    }
  }
  
  gmax00( (char *) &derror, 1, 5, 16, proclist[0], nprocs, proclist, scrat);

/*
  printf(" derror = %g normA %g ulp %g \n", derror, normA, ulp);
*/
  /*
  for ( j = 0; j < *n; j++ )
    printf(" d[%d] = %g \n", j, d[j]);
  for ( j = 0; j < *n; j++ )
    printf(" e[%d] = %g \n", j, e[j]);
  */
  
  *res = derror / normA / ulp;
  
  return;
}


