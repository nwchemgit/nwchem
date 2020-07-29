/*-----------------------------------------------------*\
 $Id$
\*-----------------------------------------------------*/
#if defined(MACX)
#include <stdio.h>
#endif
#include <stdlib.h>
#include "typesf2c.h"

void linux_sran_(Integer* input_seed)
{
  unsigned int seed;

  seed = (unsigned) *input_seed;
  (void) srandom(seed);
}
double linux_rand_(void)
{
  return (double) (((double) random())/(double) RAND_MAX);
}
