/*
 $Id: linux_cpu.c,v 1.2 1997-10-31 20:45:33 d3e129 Exp $
 */

#include <time.h>

double linux_cputime_(void)
{
  return ((double) clock()) / ((double) CLOCKS_PER_SEC);
}
