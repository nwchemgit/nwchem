/*$Id: util_date.c,v 1.4 1995-02-02 17:51:45 d3g681 Exp $*/
#include <sys/types.h>
#include <time.h>
#ifndef IPSC
#include <sys/time.h>
#endif
#ifdef CRAY
#define util_date_ UTIL_DATE
#endif
extern int string_to_fortchar(char *, int, const char *);
  
void util_date_(date)
     char *date;
/*
  Routine to return to FORTRAN the current date in 
  same format as the C routine ctime.

  
  character*(*) date
  call util_date(date)
*/
{
#ifdef CRAY
  char farray[256];
  int len,i;
#endif
  int nlen = 26;

  time_t t = time((time_t *) 0);
  char *tmp = ctime(&t);
  tmp[24] = 0;
  
#ifdef CRAY
  len = strlen(tmp);

  for (i=0; i<len; i++)
    farray[i] = tmp[i];
  for (i=len; i<nlen; i++)
    farray[i] = ' ';
  strcpy( date, tmp);
#else
  (void) string_to_fortchar(date, nlen, tmp);
#endif
}


