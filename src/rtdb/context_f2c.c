#include <stdio.h>
#include <stdlib.h>
#include <string.h>
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

static int fortchar_to_string(const char *f, int flen, char *buf, 
			      const int buflen)
{
  while (flen-- && f[flen] == ' ');
  if ((flen+1) >= buflen) return 0; /* Won't fit */
  flen++;
  buf[flen] = 0;
  while(flen--) buf[flen] = f[flen];
  return 1;
}
static int string_to_fortchar(char *f, int flen, char *buf)
{
  int len = strlen(buf), i;

  if (len > flen) return 0; /* Won't fit */
  for (i=0; i<len; i++) f[i] = buf[i];
  for (i=len; i<flen; i++) f[i] = ' ';
  return 1;
}

logical context_set_(const char *string, int len)
{
  char buf[MAX_CLEN];

  if (!fortchar_to_string(string, len, buf, sizeof buf)) {
    fprintf(stderr, "context_set: string too long? %s\n", string);
    fflush(stderr);
    return 0;
  }  
  if (context_set(buf)) {
    return FORTRAN_TRUE;
  } else {
    return FORTRAN_FALSE;
  }
}

logical context_get_(const char *string, int len)
{
  char *tmp = context_get();
  int status = string_to_fortchar(string, len, tmp);

  free(tmp);
  if (status) {
    return FORTRAN_TRUE;
  } else {
    return FORTRAN_FALSE;
  }
}

logical context_rtdb_store_(integer *prtdb)
{
  if (context_rtdb_store((int) *prtdb)) {
    return FORTRAN_TRUE;
  } else {
    return FORTRAN_FALSE;
  }
}

logical context_rtdb_load_(integer *prtdb)
{
  if (context_rtdb_load((int) *prtdb)) {
    return FORTRAN_TRUE;
  } else {
    return FORTRAN_FALSE;
  }
}
logical context_push_(const char *string, int len)
{
  char buf[MAX_CLEN];

  if (!fortchar_to_string(string, len, buf, sizeof buf)) {
    fprintf(stderr, "context_push: string too long? %s\n", string);
    fflush(stderr);
    return 0;
  }  
  if (context_push(buf)) {
    return FORTRAN_TRUE;
  } else {
    return FORTRAN_FALSE;
  }
}
logical context_pop_(const char *string, int len)
{
  char buf[MAX_CLEN];

  if (!fortchar_to_string(string, len, buf, sizeof buf)) {
    fprintf(stderr, "context_pop: string too long? %s\n", string);
    fflush(stderr);
    return 0;
  }  
  if (context_pop(buf)) {
    return FORTRAN_TRUE;
  } else {
    return FORTRAN_FALSE;
  }
}
logical context_prefix_(const char *name, char *result, int nlen, int rlen)
{
  char buf[MAX_CLEN], rbuf[MAX_CLEN];

  if (!fortchar_to_string(name, nlen, buf, sizeof(buf))) {
    fprintf(stderr,"context_rtdb_match: buffer too small\n");
    return FORTRAN_FALSE;
  }
  if (!context_prefix(buf, rbuf, sizeof rbuf)) {
    return FORTRAN_FALSE;
  }
  if (!string_to_fortchar(result, rlen, rbuf)) {
    fprintf(stderr,"context_prefix: result too small\n");
    return FORTRAN_FALSE;
  }

  return FORTRAN_TRUE;
}
logical context_rtdb_match_(integer *prtdb, const char *name, char *result,
			    int nlen, int rlen)
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
