#include <time.h>

double linux_cputime_(void)
{
  return ((double) clock()) / ((double) CLOCKS_PER_SEC);
}
