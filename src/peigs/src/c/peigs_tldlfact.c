/*
 $Id$
*/

#include <stdio.h>
#include <math.h>

#include "globalp.c.h"

void peigs_tldlfact(Integer  *n, DoublePrecision *d, DoublePrecision *e, DoublePrecision *dplus, DoublePrecision *lplus)
{
  
  Integer i,  msize=*n-1;
  
  msize = *n -1 ;
  dplus[0] = d[0];
  for (i = 0; i < msize; i++ ){
    lplus[i] = e[i+1]/dplus[i];
    dplus[i+1] = d[i+1] - lplus[i]*e[i+1];
  }
  /*  for (i = 0; i < msize; i++ ){
          printf(" %ld lp dp  %f %f \n", i, lplus[i], dplus[i]);
  }*/

  return;
}

