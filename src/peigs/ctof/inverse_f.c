/*
 $Id$
 * **************************************************
 *
 *  this is the Fortran interface to the
 *  C routine inverseL
 *
 *  matrix stored in the lower triangular
 *  packed storage format
 *
 *  on exit, it overwrites the current lower triangular matrix
 *
 */

#include <stdio.h>
#include <math.h>
#include "globalp.c.h"

/*
  
  This routine inverts a lower triangular matrix using a list-column distributed
  storage format.  On output the inverse matrix overwrites the input matrix
  with the same storage mapping.
  
  */

void inversel_( msize, map, matrix, iwork, work, info )
     Integer *msize, *map, *iwork, *info;
     DoublePrecision *matrix, *work;
     
     /*
       msize =(input) integer, size of the matrix
       
       matrix =(input) DoublePrecision precision array containing the matrix
       
       map =(input) integer array, length n, such that for the i-th column,
                      map[i] = real processor id holding this column
       
      iwork =(input) integer scratch space of size nprocs
		      
      work =(input) DoublePrecision precision array of length msize
      
       */
{
  Integer n, i, k, me;
  Integer *mapvecs, *iscrat, nvecs;
  
  DoublePrecision **columnM;
  DoublePrecision *dscrat;
  
  extern Integer mxmynd_();
  extern Integer fil_mapvec_();
  extern void xerbla_();
  extern void inverseL();
  
  me = mxmynd_();
  
  if( info == NULL ) {
    i = -6;
    xerbla_( "INVERSE ", &i);
  }
  
  if( msize == NULL ) {
    i = -1;
    xerbla_( "INVERSE ", &i);
  }
  else
    if( *msize < 1 ){
      i = -1;
      xerbla_( "INVERSE ", &i);
    }
    else {
      iscrat = map;
      for( k = 0; k < *msize; k++ ) {
	if( iscrat == NULL ) {
	  i = -3;
	  xerbla_( "INVERSE \n", &i);
	}
	else
	  iscrat++;
      }
    }
  
  /*
    at this point inputs are minimally acceptable we'll let inverse check the map distribution
    */
  
  n = *msize;  ;
  
  iscrat = iwork;
  mapvecs = iscrat;
  nvecs = fil_mapvec_( &me, &n, map, mapvecs );
  
  /*
    don't use mapvecs except in inverse
    */
  
  columnM =(DoublePrecision **) work;
  dscrat =(DoublePrecision *)( work + nvecs );
  
  k = 0;
  for( i = 0; i < nvecs; i++ ) {
    columnM[i] = &matrix[k];
    k += n - mapvecs[i];
  }
  
  inverseL( msize, columnM, map, iscrat, dscrat, info);
  return;
  
}

