#include <stdio.h>
#include <math.h>

#include "globalp.c.h"

void peigs_tldlfact(Integer  *n, DoublePrecision *d, DoublePrecision *e, DoublePrecision *dplus, DoublePrecision *lplus)
{
  
  Integer i, j, msize=*n-1;

  
  dplus[0] = d[0];
  j = 1;
  for (i = 0; i < msize; i++ ){
    lplus[i] = e[j] /dplus[i];
    dplus[j] = d[j] - lplus[i]*e[j];
    j++;
  }

  return;
}




