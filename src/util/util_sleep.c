/*
 $Id: util_sleep.c,v 1.5 2003-08-13 20:22:23 edo Exp $
 */

#if defined(CRAY) &&!defined(__crayx1)

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
