#include <sys/types.h>
#include <time.h>
#if !defined(IPSC) && !defined(WIN32)
#include <sys/time.h>
#endif

#if defined(WIN32) &&!defined(__MINGW32__)
#define util_date_ UTIL_DATE
#include "typesf2c.h"
#endif

int string_to_fortchar(char *, int, const char *);

/*
  Routine to return to FORTRAN the current date in 
  same format as the C routine ctime.

  character*(*) date
  call util_date(date)
*/


void util_date_(char *date, int nlen)
{
  time_t t = time((time_t *) 0);
  char *tmp = ctime(&t);
  
  (void) string_to_fortchar(date, nlen, tmp);
}
