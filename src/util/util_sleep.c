/*
 $Id: util_sleep.c,v 1.4 1999-11-13 03:20:15 bjohnson Exp $
 */

#if defined(CRAY)

void UTIL_SLEEP(long *t)
{}

#elif defined(WIN32)

#include "winutil.h"
#include "typesf2c.h"
void FATR UTIL_SLEEP(long *t) 
{
  unsigned s = *t;
  sleep(s);
}

#else

#include <unistd.h>
void util_sleep_(long *t) 
{
  unsigned s = *t;
  sleep(s);
}

#endif
