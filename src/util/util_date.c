#include <sys/types.h>
#include <time.h>
#ifndef IPSC
#include <sys/time.h>
#endif

extern int string_to_fortchar(char *, int, const char *);
  
void util_date_(char *date, int dlen)
/*
  Routine to return to FORTRAN the current date in 
  same format as the C routine ctime.

  
  character*(*) date
  call util_date(date)
*/
{
  time_t t = time((time_t *) 0);
  char *tmp = ctime(&t);
  tmp[24] = 0;
  
  (void) string_to_fortchar(date, dlen, tmp);
}
