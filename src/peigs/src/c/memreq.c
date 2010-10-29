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
/* ****************************************
   
   PeIGS utility memreq_ for C routines.
   
   returns the scratch memory requirements for the PDSPGV,
   PDSPEV, and PDSPTRI routines.
   
   */

#include <stdio.h>
#include <memory.h>

#include "globalp.c.h"

#define max(a,b) ((a) > (b) ? (a) : (b))

#define BDR sizeof(DoublePrecision)/sizeof(Integer)

void memreq_(type, n, mapA, mapB, mapZ, isize, rsize, ptr_size, iscratch )
     Integer *type, *n, *mapA, *mapB, *mapZ, *isize, *rsize, *ptr_size, *iscratch;
{

/*
    this subroutine computes the memory requirements for this processor to
    safely setup the eigensystem problem.  On output, isize, risize, and ptr_size
    are the sizes ( in the respective data types ) of the data needed by the eigensystem
    program.  It performs the basic error checks and assumes that the user has
    input the correct information into type, n, mapA, mapB, mapZ
    
    input:
                 memory for 

   *type = 0    the generalized eigensystem problem (pdspgvx or pdspgv )
         = 1    the standard eigensystem problem (pdspevx or pdspev )
         = 2    the tri-diagonal standard eigensystem problem (pdsptri )
         = 3    the standard eigensystem problem (no eigenvectors )
 
     n: size of the matrix
 
     mapA : integer array of size n;
     mapA[i], 0<= i < n, is the processor holding this vector
 
     mapB : integer array of size n;
     mapB[i], 0<= i < n, is the processor holding this vector
    
    mapZ : integer array of size n
    mapZ[i], 0<= i < n, is the processor holding this vector
    
    output:
    
    isize = size of integer buffer
    rsize = size of buffer for DoublePrecision precision
    ptr_size= size of double precision buffer for pointers
    
    at this point one can allocate the required memory using C
    or using off-sets from Fortran arrays
*/
  
  Integer me, msize, indx, i;
  Integer           linfo, itype, nextra;
  Integer           naproc, use_inverse, i_bufsiz, i_memreq, i_same,
                    neleA, neleA_max, neleB_max;
  Integer           i_tmp, d_tmp, d_tmp2, ppd_tmp;
  Integer           i_dstebz, i_pstebz, i_clustrf, i_inv_it,
                    i_mgs, i_clustrinv, i_pstein, i_pdsptri,
                    i_tred2, i_mxm25, i_pdspevx,
                    d_dstebz, d_pstebz, d_clustrf, d_inv_it,
                    d_mgs, d_clustrinv, d_pstein, d_pdsptri,
                    d_tred2, d_mxm25, d_pdspevx,
                    ptr_dstebz, ptr_pstebz, ptr_clustrf, ptr_inv_it,
                    ptr_mgs, ptr_clustrinv, ptr_pstein, ptr_pdsptri,
                    ptr_tred2, ptr_mxm25, ptr_pdspevx;

  Integer *iptr;
  char            msg[ 25 ];
  Integer nvecsA, nvecsQ, nvecsQ_max, nvecsZ, nvecsZ_max, nvecsB;
  

  extern Integer count_list(), ci_size_();
  extern Integer mxmynd_(), mxnprc_();
  extern void   memreq2(), l_exit_();
  extern void     mxinit_();
  
   mxinit_();

   me     = mxmynd_();
   naproc = mxnprc_();

   strcpy( msg,  "Error in memreq." );  
  
   linfo = 0;

   if ( type  == NULL )
     linfo = -1;
   else if ( n == NULL )
     linfo = -2;
   else if ( mapA == NULL )
     linfo = -3;
   else if ( mapB == NULL )
     linfo = -4;
   else if ( mapZ == NULL )
     linfo = -5;
   else if ( isize == NULL )
     linfo = -6;
   else if ( rsize == NULL )
     linfo = -7;
   else if ( ptr_size == NULL )
     linfo = -8;
   else if ( iscratch == NULL )
     linfo = -9;

   if ( linfo != 0 ) {
     l_exit_( &linfo, msg );
     return;
   }
  
   if (( *type < 0 ) ||  ( *type > 2 )) {
     linfo = -1;
     l_exit_( &linfo, msg );
     return;
   }
  
   if ( *n < 1 ) {
     linfo = -2;
     l_exit_( &linfo, msg );
     return;
   }
  
   iptr = mapA;
   for ( indx = 0; indx < *n; indx++ ) {
     if ( (iptr++ ) == NULL ) {
       linfo = -3;
       l_exit_(&linfo, msg);
       return;
     }
     if ( mapA[indx] < 0  || mapA[indx] >= naproc ) {
       linfo = -3;
       l_exit_(&linfo, msg);
       return;
     }
   }
  
   iptr = mapB;
   if ( *type == 0 ) {
     for ( indx = 0; indx < *n; indx++ ) {
       if ( (iptr++ ) == NULL ) {
         linfo = -4;
         fprintf(stderr, " me = %d error in mapB in memreq.c \n",me);
         l_exit_(&linfo,msg);
         return;
       }
       if ( mapB[indx] < 0  || mapB[indx] >= naproc ) {
         linfo = -4;
         l_exit_(&linfo, msg);
         return;
       }
     }
   }

   iptr = mapZ;
   for ( indx = 0; indx < *n; indx++ ) {
     if ( (iptr++ ) == NULL ) {
       linfo = -5;
       fprintf(stderr, " me = %d error in mapZ in memreq.c \n",me);
       l_exit_(&linfo,msg);
       return;
     }
     if ( mapZ[indx] < 0  || mapZ[indx] >= naproc ) {
       linfo = -5;
       l_exit_(&linfo, msg);
       return;
     }
   }

   itype = *type;
   msize = *n;
  
   /*
    *  nextra = some extra space for all data types, just in case
    */

   nextra = 2 * ( max( msize,  naproc ) + 20 );

   /*
    *  Use L inverse when more than 1 processor allocated.
    */

   if( naproc == 1 )
     use_inverse = 0;
   else
     use_inverse = 1;

   /*
    *  Currently code always uses L inverse.
    */

   use_inverse = 1;

   i_memreq = 0;

   /*
    *  Assume mapQ == mapA
    */

  nvecsA = count_list( me, mapA, &msize);  
  nvecsQ = nvecsA;

  if ( itype == 0 )
    nvecsB = count_list( me, mapB, &msize);  
  else
    nvecsB = 0;
  
  nvecsZ = count_list( me, mapZ, &msize);  
  
  nvecsQ_max = 1;
   for( i = 0; i < naproc; i++ ) 
     nvecsQ_max = max( nvecsQ_max, count_list( i, mapA, &msize) );

   nvecsZ_max = 1;
   for( i = 0; i < naproc; i++ ) 
     nvecsZ_max = max( nvecsZ_max, count_list( i, mapZ, &msize) );

   /* MXCOMBV1 workspace. */

   i_bufsiz = msize;
   
   /*
    *====================
    *     pdsptri
    *====================
    */

   /* dstebz */

   i_dstebz   = 10 * msize;
   d_dstebz   = 10 * msize;
   ptr_dstebz = 0;


   /* pstebz */

   i_pstebz   = msize + 9 + max( msize + naproc, i_dstebz );

   d_pstebz   = max( msize, naproc );
   d_pstebz   = max( d_pstebz, d_dstebz );
   d_pstebz   = max( msize + 7 + d_pstebz, i_bufsiz );

   ptr_pstebz = 0;


   /* clustrf */

   i_clustrf   = 2 * msize + 4;
   d_clustrf   = 0;
   ptr_clustrf = 0;


   /* inv_it */

   i_inv_it   = msize;
   d_inv_it   = 5 * msize;
   ptr_inv_it = 0;


   /* mgs */

   i_mgs   = naproc + msize + max( naproc, msize );

   d_mgs   = msize * max( 2, nvecsZ_max  + 3 );
   if( naproc == 1 ) 
     d_mgs = msize;

   ptr_mgs = 0;


   /* clustrinv */

   i_clustrinv   = max( i_mgs, i_inv_it );
   i_clustrinv   = max( i_clustrinv, naproc );

   d_clustrinv   = max( d_mgs, d_inv_it );
   d_clustrinv   = max( d_clustrinv, msize );
   if( naproc > 1 )
     d_clustrinv  += nvecsZ_max * msize;

   ptr_clustrinv = 0;


   /* pstein */

   i_pstein   = max( 2 * naproc, i_clustrf ) + 6;
   i_pstein   = 12 * msize + 7 + max( i_pstein, i_clustrinv );

   d_pstein   = max( i_bufsiz, d_clustrinv ) + 6;
   d_pstein   = msize + 1 + max( d_pstein, naproc );

   ptr_pstein = 0;


   /* pdsptri */

   i_pdsptri   = max( naproc, i_pstein );
   i_pdsptri   = 2 * msize + 16 + max( i_pdsptri, i_pstebz );

   d_pdsptri   = max( 3, d_pstein );
   d_pdsptri   = max( d_pdsptri, d_pstebz );
   d_pdsptri   = max( d_pdsptri, i_bufsiz );

   ptr_pdsptri = ptr_pstein;


   if( itype == 2 ) {
     *isize    =  i_pdsptri + nextra;
     *rsize    =  d_pdsptri + nextra;
     *ptr_size =  ptr_pdsptri + nextra;

     return;
   }

   /*
    *====================
    *     pdspevx, pdspev
    *====================
    */

   /* tred2 */

   i_tred2   = 3 * msize + 5;
   d_tred2   = msize + 1 + i_bufsiz;
   ptr_tred2 = 0;

   /*
    * mxm25
    *
    * For mxm25 assume number of eigenvalues = msize.
    */

   i_mxm25   = 3 * msize + 2 * msize;
   d_mxm25   = max( 2 * msize , ( nvecsZ + 2 * nvecsQ_max ) * msize );
   ptr_mxm25 = 0;


   /* pdspevx */

   i_pdspevx = max( naproc, i_tred2 );
   i_pdspevx = max( i_pdspevx, i_pstebz );
   i_pdspevx = max( i_pdspevx, i_pstein );
   i_pdspevx = max( i_pdspevx, i_mxm25 );
   i_pdspevx += 23 * msize; 
   i_pdspevx = max( i_pdspevx, naproc + 11 + max( msize + naproc, i_memreq ) );

   d_pdspevx = max( d_tred2, d_pstebz );
   d_pdspevx = max( d_pdspevx, d_pstein );
   d_pdspevx = max( d_pdspevx, d_mxm25 );
   d_pdspevx += ( nvecsQ + 2 )  * msize; 
   d_pdspevx += 10*msize;

   d_pdspevx = max( d_pdspevx, naproc + 7 );
   d_pdspevx = max( d_pdspevx, i_bufsiz );

   ptr_pdspevx = nvecsQ + ptr_pstein;

   /*
    * pdspev: requires no more space than pdspevx, so do nothing
    */

   if( itype == 1 ) {
      *isize    =  i_pdspevx + nextra;
      *rsize    =  d_pdspevx + nextra;
      *ptr_size =  ptr_pdspevx + nextra;

      return;
   }


   /*
    *===================================================================
    *  pdspgvx/pdspgv      (pdspgv requires no more memory than pdspgvx)
    *===================================================================
    */


   neleA = ci_size_( &me, &msize, mapA ); 

   neleA_max = 0;
   neleB_max = 0;
   for( i = 0; i < naproc; i++ )  {
     neleA_max = max( neleA_max, ci_size_( &i, &msize, mapA ) ); 
   }

  if ( itype == 0 ) {
    for( i = 0; i < naproc; i++ )  {
      neleB_max = max( neleB_max, ci_size_( &i, &msize, mapB ) ); 
    }
  }

  
  
  /* Conjugation and back-transformation of eigenvectors */
  
  if( use_inverse == 0 ) {
    
    /* Don't use L inverse */
    
    /* forwardLL */
    
    i_tmp   = 0;
    d_tmp   = msize;
    ppd_tmp = 0;
    
    /* forwardLU */
    
    i_tmp = max( i_tmp, msize );
    d_tmp = max( d_tmp, 2 * msize );
    
    /* upperUF */
    
    i_tmp = max( i_tmp, 3 * msize + naproc + 1 );
    d_tmp = max( d_tmp, msize );
    
  }
  else {
    
    /* Use L inverse */
    
    i_same = 1 ;
    for( i = 1; i < msize; i++ )
      if( mapA[i] != mapB[i] ) i_same = 0;
    
    if( i_same == 0 ) {
      
      /* Use the mapA not equal mapB path */
      
      /* pmmLUL */
      
      i_tmp   = msize + 1;
      d_tmp   = 2 * msize + neleA;
      ppd_tmp = nvecsA;
      
      /* pmmLSL */
      
      i_tmp   = max( i_tmp,   msize + 1 );
      d_tmp   = max( d_tmp,   2 * msize + neleA );
      ppd_tmp = max( ppd_tmp, nvecsA );
      
    } 
    else {
      
      
      /* Use the mapA = mapB path */
      
      /* de_sym */
      
      i_tmp   = 2 * msize;
      d_tmp   = 2 * neleA_max;
      ppd_tmp = 0;
      
      /* mxm_llx */
      
      i_tmp   = max( i_tmp,   5 * msize + 1);
      
      d_tmp2 = max( i_bufsiz, 2 * msize );
      d_tmp2 = max( d_tmp2,   nvecsA * msize + 2 * neleB_max );
      d_tmp  = max( d_tmp, d_tmp2 );
      
      ppd_tmp = max( ppd_tmp, 0 );
      
      
      /* pmmlsl2 */
      
      i_tmp   += 0;
      d_tmp   += nvecsA * msize;
      ppd_tmp += nvecsA;
      
      /* lu_mxm2 */
      
      i_tmp   = max( i_tmp, 5 * msize + 1);
      
      d_tmp2 = max( i_bufsiz, 2 * msize );
      d_tmp2 = max( d_tmp2,   2 * ( neleB_max + neleA_max ) );
      d_tmp  = max( d_tmp, d_tmp2 );
      
      ppd_tmp = max( ppd_tmp, 0 );
    }
    
    /* lsl_conjugation2 */
    
    i_tmp   += 2 * msize;
    d_tmp   += 0;
    ppd_tmp += 0;
    
    
    /* inverseL */
    
    i_tmp = max( i_tmp, 3 * msize + 3 );
    d_tmp = max( d_tmp, i_bufsiz );
    d_tmp = max( d_tmp, msize );
    
    /* mxm5x */
    
    i_tmp = max( i_tmp, 5 * msize + 1 );
    d_tmp = max( d_tmp, i_bufsiz );
    d_tmp = max( d_tmp, 2 * msize );
    d_tmp = max( d_tmp, nvecsZ * msize + 2 * neleB_max );
    
  }
  
  /* choleski */
  
  i_tmp = max( i_tmp, 3 * msize + 3 );
  d_tmp = max( d_tmp, i_bufsiz );
  d_tmp = max( d_tmp, msize + 1 );
  d_tmp += msize + msize;
  
  /* pdspgvx */
  
  i_tmp = max( i_tmp, naproc );
  i_tmp = max( i_tmp , i_pdspevx );
  i_tmp += 2 * msize + 2 * naproc + 13;
  
  i_tmp = max( i_tmp, naproc + i_memreq);
  
  d_tmp = max(  d_tmp, d_pdspevx );
  d_tmp = max(  d_tmp, i_bufsiz );
  d_tmp = max(  d_tmp, naproc + 7 );
  d_tmp += msize + msize;
  
  ppd_tmp = max(  ppd_tmp, ptr_pdspevx );

  *isize    =  i_tmp + nextra + 4*msize + 10*naproc + 16*naproc ;
  i = *isize;
  i = ( i + sizeof(DoublePrecision) ) % sizeof(DoublePrecision) ;
  *isize += *isize + i + 32*msize;
  
  
  *rsize    =  d_tmp + nextra + 8*msize + 3*(msize*msize/naproc) + 2;
  i = *rsize;
  i = ( i + sizeof(DoublePrecision) ) % sizeof(DoublePrecision) ;
  *rsize += + i + 40*msize + msize*msize/naproc;
  
  *ptr_size =  ppd_tmp + nextra + 8*msize + 20*msize;
  i = *ptr_size + 5*msize;
  *ptr_size += *ptr_size + 40*msize ;

  
  return;
}
