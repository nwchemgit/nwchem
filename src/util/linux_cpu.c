/*
 $Id: linux_cpu.c,v 1.3 1999-11-30 01:51:42 edo Exp $
 */

#include <sys/times.h>
#include <time.h>

double linux_cputime_(void)
{
  struct tms cput;

  (void) times(&cput);

  return ((clock_t) cput.tms_utime) / ((double) CLK_TCK);
}
