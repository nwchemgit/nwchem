/*$Id: context_f2c.c,v 1.4 1995-10-17 08:56:05 d3g681 Exp $*/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#ifdef CRAY
#include <fortran.h>
#endif
#include "context.h"

typedef long logical;		/* Equivalent C type to FORTRAN logical */
typedef long integer;		/* Equivalent C type to FORTRAN integer */
#define FORTRAN_TRUE  ((logical) 1)
#define FORTRAN_FALSE ((logical) 0)

#define MAX_CLEN 4096
#ifdef CRAY 
static int fortchar_to_string(_fcd f, int flen, char *buf, 
			      const int buflen)
#else
static int fortchar_to_string(const char *f, int flen, char *buf, 
			      const int buflen)
#endif
{
#ifdef CRAY
  char *fstring = _fcdtocp(f);
  flen = _fcdlen(f);

  while (flen-- && fstring[flen] == ' ')
    ;

  if ((flen+1) >= buflen)
    return 0;			/* Won't fit */

  flen++;
  buf[flen] = 0;
  while(flen--)
    buf[flen] = fstring[flen];

#else

  while (flen-- && f[flen] == ' ')
    ;

  if ((flen+1) >= buflen)
    return 0;			/* Won't fit */

  flen++;
  buf[flen] = 0;
  while(flen--)
    buf[flen] = f[flen];
#endif
  return 1;
}
static int string_to_fortchar( f, flen, buf)
  int flen;
  char *buf;
#ifdef CRAY
  _fcd f;
#else
  char *f;
#endif
{
  int len = strlen(buf), i;
#ifdef CRAY
  flen = _fcdlen(f);
#endif

  if (len > flen) 
    return 0;			/* Won't fit */

#ifdef CRAY
  for (i=0; i<len; i++)
    _fcdtocp(f)[i] = buf[i];
  for (i=len; i<flen; i++)
    _fcdtocp(f)[i] = ' ';
#else
  for (i=0; i<len; i++)
    f[i] = buf[i];
  for (i=len; i<flen; i++)
    f[i] = ' ';
#endif
  return 1;
}

#ifdef CRAY
logical context_set_(_fcd  string, int len)
#else
logical context_set_(const char *string, int len)
#endif
{
  char buf[MAX_CLEN];

  if (!fortchar_to_string(string, len, buf, sizeof buf)) {
    fprintf(stderr, "context_set: string too long? %s\n", string);
    fflush(stderr);
    return 0;
  }  
  if (context_set(buf))
    return FORTRAN_TRUE;
  else
    return FORTRAN_FALSE;
}

logical context_get_(char *string, int len)
{
  char *tmp = context_get();
  int status = string_to_fortchar(string, len, tmp);

  free(tmp);
  if (status)
    return FORTRAN_TRUE;
  else
    return FORTRAN_FALSE;
}

logical context_rtdb_store_(integer *prtdb)
{
  if (context_rtdb_store((int) *prtdb))
    return FORTRAN_TRUE;
  else
    return FORTRAN_FALSE;
}

logical context_rtdb_load_(integer *prtdb)
{
  if (context_rtdb_load((int) *prtdb))
    return FORTRAN_TRUE;
  else
    return FORTRAN_FALSE;
}
#ifdef CRAY
logical context_push_(_fcd string, int len)
#else
logical context_push_(const char *string, int len)
#endif
{
  char buf[MAX_CLEN];

  if (!fortchar_to_string(string, len, buf, sizeof buf)) {
    fprintf(stderr, "context_push: string too long? %s\n", string);
    fflush(stderr);
    return 0;
  }  
  if (context_push(buf))
    return FORTRAN_TRUE;
  else
    return FORTRAN_FALSE;
}
#ifdef CRAY
logical context_pop_(_fcd string, int len)
#else
logical context_pop_(const char *string, int len)
#endif
{
  char buf[MAX_CLEN];

  if (!fortchar_to_string(string, len, buf, sizeof buf)) {
    fprintf(stderr, "context_pop: string too long? %s\n", string);
    fflush(stderr);
    return 0;
  }  
  if (context_pop(buf))
    return FORTRAN_TRUE;
  else
    return FORTRAN_FALSE;
}
#ifdef CRAY
logical context_prefix_(_fcd name, _fcd result, int nlen, int rlen)
#else
logical context_prefix_(const char *name, char *result, int nlen, int rlen)
#endif
{
  char buf[MAX_CLEN], rbuf[MAX_CLEN];

  if (!fortchar_to_string(name, nlen, buf, sizeof(buf))) {
    fprintf(stderr,"context_rtdb_match: buffer too small\n");
    return FORTRAN_FALSE;
  }
  if (!context_prefix(buf, rbuf, sizeof rbuf))
    return FORTRAN_FALSE;
  if (!string_to_fortchar(result, rlen, rbuf)) {
    fprintf(stderr,"context_prefix: result too small\n");
    return FORTRAN_FALSE;
  }

    return FORTRAN_TRUE;
}
#ifdef CRAY  
logical context_rtdb_match_(integer *prtdb, _fcd name, _fcd result,
			    int nlen, int rlen)
#else
logical context_rtdb_match_(integer *prtdb, const char *name, char *result,
			    int nlen, int rlen)
#endif
{
  char buf[MAX_CLEN], rbuf[MAX_CLEN];

  if (!fortchar_to_string(name, nlen, buf, sizeof(buf))) {
    fprintf(stderr,"context_rtdb_match: buffer too small\n");
    return FORTRAN_FALSE;
  }
  if (!context_rtdb_match((int) *prtdb, buf, sizeof(rbuf), rbuf))
    return FORTRAN_FALSE;
  if (!string_to_fortchar(result, rlen, rbuf)) {
    fprintf(stderr,"context_rtdb_match: result too small\n");
    return FORTRAN_FALSE;
  }

  return FORTRAN_TRUE;
}
