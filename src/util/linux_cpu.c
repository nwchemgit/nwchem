/*
 $Id$
 */

#include <sys/time.h>
#include <unistd.h>
#include <stdio.h>
#if defined(__MINGW32__)
double linux_cputime_(void){
  return ((double) clock()) / ((double) CLOCKS_PER_SEC);
}
#else
#include <sys/resource.h>

double linux_cputime_(void)
{
  struct rusage rusage_out;

   (void) getrusage (RUSAGE_SELF, &rusage_out);

  return ((double)rusage_out.ru_utime.tv_usec* 1E-6 + (double)(rusage_out.ru_utime.tv_sec));
}
#endif
