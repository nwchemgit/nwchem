/*$Id: util_host.c,v 1.1 1995-12-13 01:37:44 d3g681 Exp $*/
#include <stdio.h>
#ifdef CRAY
#include <fortran.h>
#endif

extern int gethostname(char *, int);

#ifdef CRAY
extern int string_to_fortchar(_fcd, int, const char *);
void UTIL_HOSTNAME(name)
     _fcd name;
{
  int namelen = _fcdlen(name);
#else
extern int string_to_fortchar(char *, int, const char *);
void util_hostname_(name, namelen)
  char *name;
  int namelen;
{
#endif
/*
  Utility routine to return hostname to FORTRAN

  character*(*) name
  call util_hostname(name)
*/
  char buf[256];

#ifdef DELTA
  (void) string_to_fortchar(name, namelen, "delta");
#else
  if (gethostname(buf, (int) sizeof(buf)) != 0)
    buf[0] = 0;
  (void) string_to_fortchar( name, namelen, buf);
#endif
}
