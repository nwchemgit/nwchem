/*
 $Id$
*/

#include <time.h>

double linux_cputime_(void)
{
  unsigned long thetime = clock();
  double value = ((double) thetime) / ((double) CLOCKS_PER_SEC);

  /* printf("t = %ld, value = %f\n", thetime, value); */

  return value;
}
