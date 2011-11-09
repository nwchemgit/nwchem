
#include <stdio.h>
#include <math.h>

#include "globalp.c.h"

#define FABS(a) ((a) > (0.e0) ? (a) : (-a))

void r_ritz_(n, d, e, eval, map, evec, scratch, info)
     Integer *n, map[], *info;
     DoublePrecision d[], e[], eval[];
     DoublePrecision **evec;
     DoublePrecision scratch[];
{
  /*
    
    w <- Tv
    
    T = (d(0:n-1), e(0:n-2)) symmetric tridiagonal form

    assume v is of 2-norm 1

    map[i] is an integer holding the processor id holding vector
    

    
  */
  
  Integer IONE=1;
  extern Integer mxmynd_();
  DoublePrecision val, val1;
  DoublePrecision *dptr, *dscrat;
  extern DoublePrecision ddot_();
  void peigs_tri_mult();
  Integer msize, me;
  Integer i, j, k, ptr;
  
  *info = 0;
  me = mxmynd_();
  dscrat = &scratch[0];
  msize = *n;
  ptr = 0;
  
  for ( i=0; i < msize; i++ ) {
    if ( map[i] == me ) {
      dptr = evec[ptr];
      peigs_tri_mult( &msize, d, e, &evec[ptr][0], scratch);
      
      /*
	scratch = T.vec
	*/
      
      
      val1 = ddot_(&msize, scratch, &IONE, &evec[ptr][0], &IONE);
      
      /*
	val1 = rayleigh quotient
	*/
      
      
      /*
	everything should be positive definite
	*/
      
      if ( fabs(eval[i]) > 1.e-16 ) {
	/*
	  should do an absolute norm here since that is all bisection can do
	  */
	val =fabs(val1 - eval[i]);
	if ( val > 1.e-13 ){
	  printf(" PeIGS rayleigh estimates r-ritz i %d eval error val = %20.16f, eval %20.16f rel_error %20.16f \n", (int)i, val1, eval[i], val);
	  *info = -100;
	  return;
	}
      }

      /*
	local ortho check; should do a global one
	*/

      k = 0;
      for ( j = 0; j < i; j++ ){
	if ( map[j] == me ){
	  val1 = ddot_(&msize, evec[ptr], &IONE, evec[k], &IONE);
	  if ( fabs(val1) > 1.e-13 ){
	    printf(" PEIGS r-ritz ortho i %d j %d dot_val = %20.16g \n", (int)i, (int)j, val1);
	    *info = -100;
	    return;
	  }
	  k++;
	}
      }
      ptr++;
    }
  }

  /*  
  dptr11 = (Integer **) malloc(5*msize*sizeof(Integer **));
  dptr1 = (Integer *) malloc(10*msize*sizeof(Integer *));

  ortho2(n, n, evec, map, dptr11, dptr1, dscrat, &val1, info);

  free(dptr1);
  free(dptr11);

  if ( fabs(val1) > 1.e-13 ){
    printf(" in peigs check max ortho = %20.16f info = %d \n", val1, *info);
    *info = -100;
  }
  */
  
  return;
}


void peigs_tri_mult(n, d, e, v, w)
     Integer *n;
     DoublePrecision d[], e[], v[], w[];
{
  /*
    w <- Tv

    T = (d(0:n-1), e(0:n-2)) symmetric tridiagonal form

  */

  Integer i, j;
  DoublePrecision a, b, c;
  Integer msize;

  msize = *n;

  for ( i = 0; i < msize; i++ ) {
    if( msize > 1 ) {
      w[0] = d[0]*v[0]+e[0]*v[1];
      for ( j = 1; j < msize-1; j++ ){
	a = e[j-1] * v[j-1];
	b = d[j] * v[j];
	c = e[j] * v[j+1];
	w[j] = a + b + c;
      }
      w[msize-1] = d[msize-1]*v[msize-1]+e[msize-2]*v[msize-2];
    }
    else {
      w[0]  = d[0]*v[0];
    }
  }
  return;
} 



/* $Id$ */
