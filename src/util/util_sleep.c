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
