#if defined (MACX) || defined (MACX64)
  #if defined(__x86_64__)
    #include <cpuid.h>
  #endif
#else
    // sched_getcpu is non-standard and may only be found in glibc,
    // so this code is not portable to Alpine Linux, BSD,
    // or anything else that uses something else, e.g. MUSL
    #include <sched.h>
    int sched_getcpu(void);
#endif

int findmycpu_ ()
{
  int cpu = -1;

#if defined (MACX) || defined (MACX64)
  #if defined(__x86_64__)
    int cpuinfo[4];
    __cpuid_count(1, 0, cpuinfo[0], cpuinfo[1], cpuinfo[2], cpuinfo[3]);
    if ( (cpuinfo[3] & (1 << 9)) == 0 ) {
      cpu = -1;
    }
    else {
      cpu = (unsigned) cpuinfo[1] >> 24;
    }
  #endif
#else
  cpu = sched_getcpu();
#endif

  return cpu;
}



