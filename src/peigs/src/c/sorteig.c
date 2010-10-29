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


void sorteig( n, neigval, vecZ, mapZ, eval, iwork, work )
     
     Integer           *n, *neigval, mapZ[], iwork[];
     DoublePrecision   eval[], work[];
     DoublePrecision  **vecZ;

/*
 * Sort eval, vecZ mapZ so eigenalues, eval, are in increasing
 * order
 *
 *    n ........  length of eigenvectors (dimension of matrices)
 *
 *    neigval ... number of eigenvalues and eigenvectors found and
 *                to be sorted
 *
 *    iwork ..... work array of length neigval + max( neigval, mxnprc_()).
 *
 *    work ...... work array of length max( n, neigval )
 */

{
    static Integer    IONE = 1;

    Integer           me, naproc, nvecsZ, indx, jj, isaved, jlast, itmp,
                      next, ileft, iright, meval, k;
    Integer           *iorder, *i_work;

    extern void      dshellsort2_(), dcopy_();
    extern Integer   qqsort();
    extern Integer   mxmynd_(), mxnprc_();


    me     = mxmynd_();
    naproc = mxnprc_();

    meval = *neigval;

    iorder = iwork;
    i_work = iwork + meval;


    /*
     * Sort eigenvalues in increasing order, but keep track of the
     * change in indices.  iorder[indx] holds the old eigenvalue index
     * and determine the permutation which mapZ and vecZ must undergo.
     */

    for ( indx = 0; indx < meval; indx++ )
       iorder[indx] = indx;

    dshellsort2_( neigval, eval, iorder );

    /* Sort equal eigenvalues by block.
     * Assumes original eigenvalues were sorted by block.
     */

    next = 0;
    ileft = next;
    iright = ileft;
    for ( indx = 0; indx < meval; indx++ ){
      if ( eval[next] == eval[ileft] ){
	next++;
      }
      else {
	iright = next-1;
	if (ileft != iright ) {
	  qqsort(iorder, ileft, iright);
	}
	ileft = next;
	next++;
      }
      
      if ( next == meval ) {
	iright = next-1;
	if ( iright != ileft ) {
	  qqsort(iorder, ileft, iright);
	}
      }
    }

    if (NO_EVEC)
      goto END;

    
    /*
     * Compute required permutation of vecZ.
     * Actually only do part of work now, i.e., part requiring
     * original mapZ. The rest is done after permuting mapZ,
     * this saves some workspace. 
     *
     * Let indx be the global index of an eigenvalue that this
     * processor owns, then set i_work[ indx ] equal to the index
     * in vecZ of the eigenvector associated with eigenvalue index.
     */

    for ( indx = 0; indx < meval; indx++ )
      i_work[ indx ] = -1;

    nvecsZ = 0;
    for ( indx = 0; indx < meval; indx++ )
      if( mapZ[indx] == me ) {
        i_work[ indx ] =  nvecsZ;
        nvecsZ++;
      }

    /*
     * Rearrange mapZ to agree with sorted eigenvalues 
     * using minimal workspace and data movement
     */

    for ( indx = 0; indx < meval; indx++ ) {

       jj = iorder[ indx ];

       if( jj >= 0  && jj != indx ) {

          isaved = indx;
          jlast  = indx;
          itmp   = mapZ[ indx ];
           
          while ( jj != isaved ) {

             jj = iorder[ jlast ];

             iorder[ jlast ] = -(jj+1);

             if( jj == isaved )
               mapZ[ jlast ] = itmp;
             else
               mapZ[ jlast ] = mapZ[ jj ];
  
             jlast = jj;
          }
       }
    }

    /*
     * Restore iorder.
     */

    for ( indx = 0; indx < meval; indx++ )
      if( iorder[ indx ] < 0 ) 
        iorder[ indx ] = -iorder[ indx ] - 1;

    /*
     * Finish computing permutation required of vecZ
     * by overwriting mapZ iorder with vecZ iorder
     */

    k = 0;
    for ( indx = 0; indx < meval; indx++ ) 
       if( i_work[ iorder[indx] ] > -1 ) {
         iorder[ k ] = i_work[ iorder[indx] ];
         k++;
       }

    /*
     * Permute vecZ to agree with sorted eigenvalues 
     * using minimal workspace and data movement.
     */

    for ( indx = 0; indx < nvecsZ; indx++ ) {

       jj = iorder[ indx ];

       if( jj >= 0  && jj != indx ) {

          isaved = indx;
          jlast  = indx;
          dcopy_( n, vecZ[ indx ], &IONE, work, &IONE );
           
          while ( jj != isaved ) {

             jj = iorder[ jlast ];

             iorder[ jlast ] = -1;

             if( jj == isaved )
               dcopy_( n, work, &IONE, vecZ[ jlast ], &IONE );
             else
               dcopy_( n, vecZ[ jj ], &IONE, vecZ[ jlast ], &IONE );
  
             jlast = jj;
          }
       }
    }

END:
    
    return;
}
    
