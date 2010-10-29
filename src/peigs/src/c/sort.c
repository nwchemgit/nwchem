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

#define min(a,b) ( (a) < (b) ? (a) : (b))

static Integer G_count( m, n, array )
     Integer m, n, *array;
{
  Integer i, j, *k;
  
  k = array;
  i = 0;
  
  for ( j = 0; j < n; j++)
    if ( array[j] == m ) i++;
  
  return(i);
}

static Integer G_index(ishift, i )
     Integer *ishift, i;
{
  Integer junk, j;
  
  junk = 0;
  for ( j = 0; j < i ; j++ )
    junk += ishift[j];
  return(junk);
}

void sort_(m, n, nsplit, isplit, iblock, iwork, eval, work)
     Integer *m, *n, *nsplit, *isplit, *iblock, *iwork;
     DoublePrecision *eval, *work;
{
  /*
    index sort
    
    sort a list of n eigenvalues according to a corresponding block number:
    e.g.: e-vals in block number 1 get sorted from small to large ...
    m is the dimension of the matrix from which the eigenvalues came.
    
    */
  
  Integer i, indx, indx2, indx3, indx4;
  DoublePrecision *dummy;
  Integer *endptr, indx5;
  extern void heapsort_ ();
  
  /*
    special case: all the blocks are of size 1 and *n = *m.
    */
  
  if ( *n == *nsplit && *m == *n) {
    for ( i = 0; i < *n; i++ ){
      indx = iblock[i]-1;
      iwork[indx] = iblock[i];
      work[indx] = eval[i];
    }
    
    for ( i = 0; i < *n; i++ ){
      iblock[i] = iwork[i];
      eval[i]   = work[i];
    }
  }
  else
    {
      
      /*
	general case
	0 < iwork < nsplit are thes temporary array for the sizes of each block
	iwork contains the sizes of each of the blocks
	*/
      
      for ( i = 0; i < *nsplit ; i++ ){
	iwork[i] = G_count(i+1, *n, &iblock[0]);
      }
      
      dummy = (DoublePrecision *) ( work + *nsplit + 1);
      
      indx2 = -1;
      for ( i = 0; i < *nsplit; i++)
	for ( indx = 0; indx < iwork[i] ; indx++ ) {
	  indx2++;
	  dummy[ indx2 ] = (DoublePrecision) 0.0e0;
	}
      
      /*
	for ( i = 0; i < *nsplit; i++ )
	dummy[i] = (DoublePrecision *) malloc( iwork[i] * sizeof(DoublePrecision));
	
	endptr = (Integer *) malloc( (*nsplit + 1)* sizeof(Integer));
	*/
      
      
      endptr = (Integer *) work;
      
      for ( i = 0; i < *nsplit+1; i++ )
	endptr[i] = -1;
      
      for ( i = 0; i < *n; i++ ){
	indx2 = iblock[i]-1;
	++endptr[indx2];
	indx = endptr[indx2];
	indx3 = G_index( iwork, indx2 );
	dummy[indx3 + indx] = eval[i];
      }
      
      /*
	fprintf(stderr, " got here iwork %d n = %d \n", iwork[0], *n);
	*/
      
      
      /*
	on the i860 dstebz doesn't always return the correct block information 
	we sort dummy[indx][*] for this case
	*/
      
      /*
	fprintf(stderr, " just before heapsort \n " );
	*/
      
      for ( i = 0; i < *nsplit ; i++ )  {
	indx3 = G_index( iwork, i );
	heapsort_ (&iwork[i], &dummy[indx3]);
      }
      
      indx3 = 0;
      for ( indx = 0; indx < *nsplit; indx++ )
	{
	  indx2 = indx+1;
	  for ( i = 0; i < iwork[indx]; i++){
	    indx4 = indx3+i;
	    iblock[indx4] = indx2;
	    indx5 = G_index( iwork, indx );
	    eval[indx4] = dummy[indx5 + i];
	  }
	  indx3 = indx3 + iwork[indx];
	}
    }
  return;
}
