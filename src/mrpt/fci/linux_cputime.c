/*
 $Id: linux_cputime.c,v 1.2 1999-07-28 00:23:43 d3e129 Exp $
*/

#include <time.h>

double linux_cputime_(void)
{
  unsigned long thetime = clock();
  double value = ((double) thetime) / ((double) CLOCKS_PER_SEC);

  /* printf("t = %ld, value = %f\n", thetime, value); */

  return value;
}
