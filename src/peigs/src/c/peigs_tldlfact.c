/*
 $Id: peigs_tldlfact.c,v 1.6 1999-11-04 22:35:26 d3g270 Exp $
*/

#include <stdio.h>
#include <math.h>

#include "globalp.c.h"

void peigs_tldlfact(Integer  *n, DoublePrecision *d, DoublePrecision *e, DoublePrecision *dplus, DoublePrecision *lplus)
{
  
  Integer i, j, msize=*n-1;
  
  msize = *n -1 ;
  dplus[0] = d[0];
  j = 1;
  for (i = 0; i < msize; i++ ){
    lplus[i] = e[j]/dplus[i];
    dplus[j] = d[j] - lplus[i]*e[j];
    j++;
  }

  return;
}

