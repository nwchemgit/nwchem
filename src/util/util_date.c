/*$Id$*/
#include <sys/types.h>
#include <time.h>
#if !defined(IPSC) && !defined(WIN32)
#include <sys/time.h>
#endif

#if defined(CRAY) && !defined(__crayx1)
#define util_date_ UTIL_DATE
#include <fortran.h>
#define FATR
#endif
#if defined(WIN32) &&!defined(__MINGW32__)
#define util_date_ UTIL_DATE
#include "typesf2c.h"
#endif

#if defined(USE_FCD)
int string_to_fortchar(_fcd, int, const char *);
#else
int string_to_fortchar(char *, int, const char *);
#endif

/*
  Routine to return to FORTRAN the current date in 
  same format as the C routine ctime.

  character*(*) date
  call util_date(date)
*/


#if defined(USE_FCD)
void FATR util_date_(_fcd date)
{
  int  nlen = _fcdlen(date);
#else
void util_date_(char *date, int nlen)
{
#endif
  time_t t = time((time_t *) 0);
  char *tmp = ctime(&t);
  
  (void) string_to_fortchar(date, nlen, tmp);
}
