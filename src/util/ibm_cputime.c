/*$Id*/
#include <time.h>

double ibm_cputime_()
{
  return clock()/CLOCKS_PER_SEC;
}
