/*
 $Id: util_sleep.c,v 1.3 1997-10-31 20:45:39 d3e129 Exp $
 */

#ifdef CRAY
void UTIL_SLEEP(long *t)
{}
#else
#include <unistd.h>
void util_sleep_(long *t) 
{
  unsigned s = *t;
  sleep(s);
}
#endif
