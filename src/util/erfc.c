#include <math.h>

extern double erfc (double);

double derfc_ (double * x)
{
  return (erfc (*x));
}

double erfc_ (double * x)
{
  return (erfc (*x));
}
/* $Id$ */
