/*
 $Id: win32_cpu.c,v 1.1 2000-02-08 23:08:14 bjohnson Exp $
*/

#include <time.h>
#include "typesf2c.h"

double FATR WIN32_CPUTIME(void)
{
  return ((double) clock()) / ((double) CLOCKS_PER_SEC);
}
