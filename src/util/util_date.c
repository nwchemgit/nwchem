/*$Id: util_date.c,v 1.6 1999-11-13 03:16:46 bjohnson Exp $*/
#include <sys/types.h>
#include <time.h>
#if !defined(IPSC) && !defined(WIN32)
#include <sys/time.h>
#endif

#ifdef CRAY
#define util_date_ UTIL_DATE
#include <fortran.h>
#endif
#ifdef WIN32
#define util_date_ UTIL_DATE
#include "typesf2c.h"
#endif

#if defined(CRAY) || defined(USE_FCD)
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


#if defined(CRAY) || defined(USE_FCD)
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
