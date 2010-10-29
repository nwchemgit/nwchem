/*-----------------------------------------------------*\
 $Id$
\*-----------------------------------------------------*/
#include <stdlib.h>
void linux_sran_(int *input_seed)
{
  unsigned int seed;

  seed = (unsigned) *input_seed;
  (void) srandom(seed);
}
double linux_rand_(void)
{
  return (double) (((double) random())/(double) RAND_MAX);
}
