/*$Id: rtdb_f2c.c,v 1.8 1995-02-02 23:22:11 d3g681 Exp $*/
#include <stdio.h>
#include <string.h>
#include "rtdb.h"
#include "macdecls.h"
#ifdef CRAY
#include "fortran.h"
#endif

typedef long logical;		/* Equivalent C type to FORTRAN logical */

typedef long integer;		/* Equivalent C type to FORTRAN integer */

#define FORTRAN_TRUE  ((logical) 1)
#define FORTRAN_FALSE ((logical) 0)


int fortchar_to_string( f,  flen, buf, buflen)
    int flen;
    char *buf;
    const int buflen;
#ifdef CRAY
    _fcd	f;		/* FORTRAN character descriptor */
#else
    const char *f;
#endif
{
#ifdef CRAY
  char *fstring;
    fstring = _fcdtocp(f);
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

int string_to_fortchar( f, flen, buf)
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


logical rtdb_parallel_(const logical *mode)
{
  /* This causes problems on machines where true != 1 (i.e. intel)
   * so it is better just to pass what we are given
   *
   * int new = (*mode == FORTRAN_TRUE); 
   */
  int new = *mode;
  int old = rtdb_parallel(new);

  if (old)
    return FORTRAN_TRUE;
  else
    return FORTRAN_FALSE;
}

#ifdef CRAY
logical rtdb_open_(_fcd filename, _fcd mode, integer *handle,
		   const integer flen, const integer mlen)
#else
logical rtdb_open_(const char *filename, const char *mode, integer *handle,
		   const integer flen, const integer mlen)
#endif
{
  char fbuf[256], mbuf[256];
  int hbuf;

  if (!fortchar_to_string(filename, flen, fbuf, sizeof(fbuf))) {
    (void) fprintf(stderr, "rtdb_open: fbuf is too small, need=%d\n",
		   (int) flen);
    return FORTRAN_FALSE;
  }

  if (!fortchar_to_string(mode, mlen, mbuf, sizeof(mbuf))) {
    (void) fprintf(stderr, "rtdb_open: mbuf is too small, need=%d\n",
		   (int) mlen);
    return FORTRAN_FALSE;
  }

  if (rtdb_open(fbuf, mbuf, &hbuf)) {
    *handle = (integer) hbuf;
    return FORTRAN_TRUE;
  }
  else {
    return FORTRAN_FALSE;
  }
}
#ifdef CRAY
logical rtdb_close_(const integer *handle, _fcd mode, const int mlen)
#else
logical rtdb_close_(const integer *handle, const char *mode, const int mlen)
#endif
{
  char mbuf[256];
  int hbuf = (int) *handle;

  if (!fortchar_to_string(mode, mlen, mbuf, sizeof(mbuf))) {
    (void) fprintf(stderr, "rtdb_close: mbuf is too small, need=%d\n", mlen);
    return FORTRAN_FALSE;
  }

  if (rtdb_close(hbuf, mbuf))
    return FORTRAN_TRUE;
  else
    return FORTRAN_FALSE;
}
#ifdef CRAY
logical rtdb_get_info_(const integer *handle, _fcd name, 
		       integer *ma_type, integer *nelem, _fcd date,
		       const int nlen, const int dlen)
#else
logical rtdb_get_info_(const integer *handle, const char *name, 
		       integer *ma_type, integer *nelem, char *date,
		       const int nlen, const int dlen)
#endif
{
  int hbuf = (int) *handle;
  char dbuf[26], nbuf[256];
  int nelbuf, typebuf;

  if (!fortchar_to_string(name, nlen, nbuf, sizeof(nbuf))) {
    (void) fprintf(stderr, "rtdb_get_info: nbuf is too small, need=%d\n", 
		   nlen);
    return FORTRAN_FALSE;
  }

  if (dlen < 24) {
    (void) fprintf(stderr, "rtdb_get_info: date must be > character*24\n");
    return FORTRAN_FALSE;
  }
    
  if (rtdb_get_info(hbuf, nbuf, &typebuf, &nelbuf, dbuf)) {
    *ma_type = (integer) typebuf;
    *nelem   = (integer) nelbuf;

    if (typebuf == MT_CHAR)	/* Fortran is ignorant of trailing null char */
      *nelem = *nelem - 1;

    if (!string_to_fortchar(date, dlen, dbuf)) {
      (void) fprintf(stderr, "rtdb_get_info: nbuf is too small, need=%d\n", 
		     nlen);
      return FORTRAN_FALSE;
    }

    return FORTRAN_TRUE;
  }
  else {
    return FORTRAN_FALSE;
  }
}
#ifdef CRAY
logical rtdb_put_(const integer *handle, _fcd name, const integer *ma_type,
		  const integer *nelem, const void *array, const int nlen)
#else
logical rtdb_put_(const integer *handle, const char *name, const integer *ma_type,
		  const integer *nelem, const void *array, const int nlen)
#endif
{
  int hbuf = (int) *handle;
  char nbuf[256];
  int nelbuf;
  int typebuf;

  if (!fortchar_to_string(name, nlen, nbuf, sizeof(nbuf))) {
    (void) fprintf(stderr, "rtdb_put: nbuf is too small, need=%d\n", 
		   nlen);
    return FORTRAN_FALSE;
  }

  nelbuf = (int) *nelem;
  typebuf= (int) *ma_type;

#ifdef DEBUG
  printf("put: rtdb=%d, mat=%d, nel=%d, name=%s\n", hbuf, typebuf, nelbuf, nbuf);
  fflush(stdout);
#endif

  if (rtdb_put(hbuf, nbuf, typebuf, nelbuf, array))
    return FORTRAN_TRUE;
  else
    return FORTRAN_FALSE;
}
#ifdef CRAY
logical rtdb_get_(const integer *handle, _fcd name, 
		  const integer *ma_type, const integer *nelem, 
		  void *array, const int nlen)
#else
logical rtdb_get_(const integer *handle, const char *name, 
		  const integer *ma_type, const integer *nelem, 
		  void *array, const int nlen)
#endif
{
  int hbuf = (int) *handle;
  char nbuf[256];
  int nelbuf;
  int typebuf;

  if (!fortchar_to_string(name, nlen, nbuf, sizeof(nbuf))) {
    (void) fprintf(stderr, "rtdb_get: nbuf is too small, need=%d\n", 
		   nlen);
    return FORTRAN_FALSE;
  }

  nelbuf = (int) *nelem;
  typebuf= (int) *ma_type;

#ifdef DEBUG
  printf("get: rtdb=%d, mat=%d, nel=%d, name=%s\n", hbuf, typebuf, nelbuf, nbuf);
  fflush(stdout);
#endif

  if (rtdb_get(hbuf, nbuf, typebuf, nelbuf, array)) {
    return FORTRAN_TRUE;
  }
  else
    return FORTRAN_FALSE;
}
#ifdef CRAY
logical rtdb_ma_get_(const integer *handle, _fcd name, integer *ma_type,
		     integer *nelem, integer *ma_handle, const int nlen)
#else
logical rtdb_ma_get_(const integer *handle, const char *name, integer *ma_type,
		     integer *nelem, integer *ma_handle, const int nlen)
#endif
{
  int hbuf = (int) *handle;
  char nbuf[256];
  int nelbuf;
  int typebuf;
  int handbuf;

  if (!fortchar_to_string(name, nlen, nbuf, sizeof(nbuf))) {
    (void) fprintf(stderr, "rtdb_ma_get: nbuf is too small, need=%d\n", 
		   nlen);
    return FORTRAN_FALSE;
  }

  if (rtdb_ma_get(hbuf, nbuf, &typebuf, &nelbuf, &handbuf)) {
    *ma_type   = (integer) typebuf;
    *ma_handle = (integer) handbuf;
    *nelem     = (integer) nelbuf;

    return FORTRAN_TRUE;
  }
  else
    return FORTRAN_FALSE;
}

logical rtdb_print_(const integer *handle, const logical *print_values)
{
  int hbuf = (int) *handle;
  int pbuf = (*print_values == FORTRAN_TRUE);

  if (rtdb_print(hbuf, pbuf))
    return FORTRAN_TRUE;
  else
    return FORTRAN_FALSE;
}
#ifdef CRAY
logical rtdb_cput_(const integer *handle, _fcd name,
		   const integer *nelem,
		   _fcd array, const int nlen, int alen)
#else
logical rtdb_cput_(const integer *handle, const char *name,
		   const integer *nelem,
		   const char *array, const int nlen, const int alen)
#endif
/*
  Insert an array of Fortran character variables into the data base.
  Each array element is striped of trailing blanks, terminated with CR,
  and appended to the list. The entire array must fit into abuf.
*/
{
  int hbuf = (int) *handle;
  char nbuf[256];
  char abuf[10240];
  int nelbuf;
  int typebuf;
  int i, left;
  char *next;
#ifdef CRAY
    alen= 0;
  for (i=0, left=sizeof(abuf), next=abuf;
       i<*nelem;
       i++, _fcdtocp(array)+=alen) {
#else
  for (i=0, left=sizeof(abuf), next=abuf;
       i<*nelem;
       i++, array+=alen) {
#endif
    
    if (!fortchar_to_string(array, alen, next, left)) {
      (void) fprintf(stderr, "rtdb_cput: abuf is too small, need=%d\n", 
		     (int) (alen + sizeof(abuf) - left));
      return FORTRAN_FALSE;
    }
#ifdef CRAY
    alen= _fcdlen(array);
#endif
    left -= strlen(next) + 1;
    next += strlen(next) + 1;
    if (i != (*nelem - 1))
      *(next-1) = '\n';
  }

  if (!fortchar_to_string(name, nlen, nbuf, sizeof(nbuf))) {
    (void) fprintf(stderr, "rtdb_cput: nbuf is too small, need=%d\n", 
		   nlen);
    return FORTRAN_FALSE;
  }

  nelbuf = strlen(abuf) + 1;
  typebuf= (int) MT_CHAR;

#ifdef DEBUG
  printf("cput: rtdb=%d, mat=%d, nel=%d, name=%s\n", hbuf, typebuf, nelbuf, nbuf);
  fflush(stdout);
#endif

  if (rtdb_put(hbuf, nbuf, typebuf, nelbuf, abuf))
    return FORTRAN_TRUE;
  else
    return FORTRAN_FALSE;
}
#ifdef CRAY
logical rtdb_cget_(const integer *handle, _fcd name,
		   const integer *nelem,
		   _fcd array, const int nlen,  int alen)
#else
logical rtdb_cget_(const integer *handle, char const *name,
		   const integer *nelem,
		   char *array, const int nlen, const int alen)
#endif
{
  int hbuf = (int) *handle;
  char nbuf[256];
  char abuf[10240];
  int nelbuf;
  int typebuf;
  int i;
  char *next;
#ifdef CRAY
  int flen=_fcdlen(array);
#endif

  if (!fortchar_to_string(name, nlen, nbuf, sizeof(nbuf))) {
    (void) fprintf(stderr, "rtdb_cget: nbuf is too small, need=%d\n", 
		   nlen);
    return FORTRAN_FALSE;
  }

  nelbuf = sizeof(abuf);
  typebuf= (int) MT_CHAR;

  if (!rtdb_get(hbuf, nbuf, typebuf, nelbuf, abuf))
    return FORTRAN_FALSE;
#ifdef CRAY
  alen =0;

  for (i=0, next=strtok(abuf, "\n");
       next;
       i++, _fcdtocp(array)+=alen, next=strtok((char *) 0, "\n")) {
#else
  for (i=0, next=strtok(abuf, "\n");
       next;
       i++, array+=alen, next=strtok((char *) 0, "\n")) {
#endif
    if (i == *nelem) {
      (void) fprintf(stderr, "rtdb_cget: array has too few elements\n");
      return FORTRAN_FALSE;
    }

    if (!string_to_fortchar(array, alen, next)) {
      (void) fprintf(stderr, "rtdb_cget: array element is too small\n");
      return FORTRAN_FALSE;
    }
#ifdef CRAY
      alen = flen;
#endif
  }
  return FORTRAN_TRUE;
}

logical rtdb_first_(const integer *handle, char *name, int nlen)
{
  char nbuf[256];
  
  if (rtdb_first((int) *handle, (int) sizeof(nbuf), nbuf) &&
      string_to_fortchar(name, nlen, nbuf))
    return FORTRAN_TRUE;
  else
    return FORTRAN_FALSE;
}

logical rtdb_next_(const integer *handle, char *name, int nlen)
{
  char nbuf[256];
  
  if (rtdb_next((int) *handle, (int) sizeof(nbuf), nbuf) &&
      string_to_fortchar(name, nlen, nbuf))
    return FORTRAN_TRUE;
  else
    return FORTRAN_FALSE;
}
#ifdef CRAY
logical rtdb_delete_(const integer *handle, _fcd name, const int nlen)
#else
logical rtdb_delete_(const integer *handle, const char *name, const int nlen)
#endif
{
  int hbuf = (int) *handle;
  char nbuf[256];

  if (!fortchar_to_string(name, nlen, nbuf, sizeof(nbuf))) {
    (void) fprintf(stderr, "rtdb_delete: nbuf is too small, need=%d\n", 
		   nlen);
    return FORTRAN_FALSE;
  }

  if (rtdb_delete(hbuf, nbuf))
    return FORTRAN_TRUE;
  else
    return FORTRAN_FALSE;
}

