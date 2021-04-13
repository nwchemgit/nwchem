#include <utmpx.h>

int sched_getcpu();

int findmycpu_ ()
{
  int cpu;
  cpu = sched_getcpu();
  return cpu;
}
