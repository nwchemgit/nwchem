/*$Id$*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#if defined(CRAY) && !defined(__crayx1)
#include <fortran.h>
#endif
#if defined(WIN32) &&!defined(__MINGW32__)
#include "typesf2c.h"
#endif
#include "context.h"

typedef long logical;		/* Equivalent C type to FORTRAN logical */
#ifdef EXT_INT
#ifdef WIN64
typedef long long integer;		/* Equivalent C type to FORTRAN integer */
#else
typedef long integer;		/* Equivalent C type to FORTRAN integer */
#endif
#else
typedef int integer;		/* Equivalent C type to FORTRAN integer */
#endif
#define FORTRAN_TRUE  ((logical) 1)
#define FORTRAN_FALSE ((logical) 0)

#ifndef WIN32
#define FATR 
#endif

#define MAX_CLEN 4096
#if (defined(CRAY) || defined(USE_FCD)) && !defined(__crayx1)
static int fortchar_to_string(_fcd f, int flen, char *buf, 
			      const int buflen)
#else
static int fortchar_to_string(const char *f, int flen, char *buf, 
			      const int buflen)
#endif
{
#if (defined(CRAY) || defined(USE_FCD))&& !defined(__crayx1)
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
#if (defined(CRAY) || defined(USE_FCD)) && !defined(__crayx1)
  _fcd f;
#else
  char *f;
#endif
{
  int len = strlen(buf), i;
#if (defined(CRAY) || defined(USE_FCD)) && !defined(__crayx1)
  flen = _fcdlen(f);
#endif

  if (len > flen) 
    return 0;			/* Won't fit */

#if (defined(CRAY) || defined(USE_FCD)) && !defined(__crayx1)
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

#if (defined(CRAY) || defined(USE_FCD)) && !defined(__crayx1)
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

#if (defined(CRAY) || defined(USE_FCD)) && !defined(__crayx1)
logical FATR context_get_(_fcd string, int len)
#else
logical context_get_(const char *string, int len)
#endif
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
#if (defined(CRAY) || defined(USE_FCD)) && !defined(__crayx1)
logical FATR context_push_(_fcd string)
#else
logical context_push_(const char *string, int len)
#endif
{
#if (defined(CRAY) || defined(USE_FCD)) && !defined(__crayx1)
  int len = _fcdlen(string);
#endif
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
#if (defined(CRAY) || defined(USE_FCD)) && !defined(__crayx1)
logical FATR context_pop_(_fcd string)
#else
logical context_pop_(const char *string, int len)
#endif
{
#if (defined(CRAY) || defined(USE_FCD)) && !defined(__crayx1)
  int len = _fcdlen(string);
#endif
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
#if (defined(CRAY) || defined(USE_FCD)) && !defined(__crayx1)
logical FATR context_prefix_(_fcd name, _fcd result)
#else
logical context_prefix_(const char *name, char *result, int nlen, int rlen)
#endif
{
#if (defined(CRAY) || defined(USE_FCD)) && !defined(__crayx1)
  int nlen = _fcdlen(name);
  int rlen = _fcdlen(result);
#endif
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
#if (defined(CRAY) || defined(USE_FCD)) && !defined(__crayx1)
logical FATR context_rtdb_match_(integer *prtdb, _fcd name, _fcd result)
#else
logical context_rtdb_match_(integer *prtdb, const char *name, char *result,
			    int nlen, int rlen)
#endif
{
#if (defined(CRAY) || defined(USE_FCD)) && !defined(__crayx1)
  int nlen = _fcdlen(name);
  int rlen = _fcdlen(result);
#endif
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
