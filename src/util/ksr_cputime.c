/* $Id: ksr_cputime.c,v 1.1 1995-07-10 17:02:12 d3e129 Exp $ */

#define USE_GETR
#if defined(USE_GETR)
/*-----------------------------*\
    use getrusage function
\*-----------------------------*/
#define __KSR_CPUTIME__
#include <stdio.h>
#include <sys/types.h>
#include <sys/time.h>
#include <sys/resource.h>

double ksr_cputime_()
{
  struct rusage myrusage;
  int errcode;
  double time_in_seconds;
  
  errcode = getrusage(RUSAGE_SELF, &myrusage);
  if (errcode != 0)
    {
      (void) fprintf(stderr, " getrusage error \n");
      return (double) -565.6589;
    }
/*
#  (void) printf(" user sec  -> %d \n",myrusage.ru_utime.tv_sec);
#  (void) printf(" user usec -> %d \n",myrusage.ru_utime.tv_usec);
#  (void) printf(" sys  sec  -> %d \n",myrusage.ru_stime.tv_sec);
#  (void) printf(" sys  usec -> %d \n",myrusage.ru_stime.tv_usec);
*/
  time_in_seconds =  (double) myrusage.ru_utime.tv_sec;
  time_in_seconds += ((double) myrusage.ru_utime.tv_usec)/(double)1000000.0;
/*debug:  (void) printf("time_in_seconds(1) = %f \n",time_in_seconds);*/
  time_in_seconds += (double) myrusage.ru_stime.tv_sec;
  time_in_seconds += ((double) myrusage.ru_stime.tv_usec)/(double)1000000.0;
/*debug:  (void) printf("time_in_seconds(2) = %f \n",time_in_seconds);*/
  return time_in_seconds;
}
#endif

#if defined(USE_CLOCK)
#define __KSR_CPUTIME__
#include <stdio.h>
#include <time.h>
double ksr_cputime_()
{
  double value;

  value = (double)clock();

/*neglect-start*\
  (void) printf("\nraw value: %f     ",(double) value);
\*neglect-end*/

  value /= CLOCKS_PER_SEC;

/*neglect-start*\
  (void) printf("sec_value: %f\n\n",(double) value);
  (void) fflush(stdout);
\*neglect-end*/
  return (double) value;
}
#endif


#ifndef __KSR_CPUTIME__
#include <stdio.h>
double ksr_cputime_()
{
  return (double) 0.0;
}
#endif
