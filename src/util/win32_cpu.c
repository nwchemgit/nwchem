/*
 $Id: win32_cpu.c,v 1.2 2001-06-28 00:16:13 edo Exp $
*/

#include <time.h>
#include "typesf2c.h"
#ifdef WIN32
double FATR WIN32_CPUTIME(void)
#else
double win32_cputime_(void)
#endif
{
  return ((double) clock()) / ((double) CLOCKS_PER_SEC);
}
