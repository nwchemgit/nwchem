#if defined (MACX) || defined (MACX64)
#include <cpuid.h>
#include <stdint.h>
#else
#include <utmpx.h>
int sched_getcpu();
#endif


int findmycpu_ ()
{
  int cpu;
#if defined (MACX) || defined (MACX64)
  int cpuinfo[4];
  __cpuid_count(1, 0, cpuinfo[0], cpuinfo[1], cpuinfo[2], cpuinfo[3]);
  if ( (cpuinfo[3] & (1 << 9)) == 0 ) {
    cpu = -1;
  }
  else {
    cpu = (unsigned) cpuinfo[1] >> 24;
  }
#else
  cpu = sched_getcpu();
#endif

  return cpu;
}



