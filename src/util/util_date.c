/*$Id: util_date.c,v 1.5 1997-02-27 00:15:52 d3g681 Exp $*/
#include <sys/types.h>
#include <time.h>
#ifndef IPSC
#include <sys/time.h>
#endif

#ifdef CRAY
#define util_date_ UTIL_DATE
#include <fortran.h>
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


#ifdef CRAY
void util_date_(_fcd date)
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
