#include <stdio.h>
#include <string.h>
#include "rtdb.h"
#include "macdecls.h"

typedef long logical;		/* Equivalent C type to FORTRAN logical */

typedef long integer;		/* Equivalent C type to FORTRAN integer */

#define FORTRAN_TRUE  ((logical) 1)
#define FORTRAN_FALSE ((logical) 0)


static int fortchar_to_string(const char *f, int flen, char *buf, 
			      const int buflen)
{
  while (flen-- && f[flen] == ' ')
    ;

  if ((flen+1) >= buflen)
    return 0;			/* Won't fit */

  flen++;
  buf[flen] = 0;
  while(flen--)
    buf[flen] = f[flen];

  return 1;
}

static int string_to_fortchar(char *f, int flen, const char *buf)
{
  int len = strlen(buf), i;

  if (len > flen) 
    return 0;			/* Won't fit */

  for (i=0; i<len; i++)
    f[i] = buf[i];
  for (i=len; i<flen; i++)
    f[i] = ' ';

  return 1;
}


logical rtdb_parallel_(const logical *mode)
{
  int new = (*mode == FORTRAN_TRUE);
  int old = rtdb_parallel(new);

  if (old)
    return FORTRAN_TRUE;
  else
    return FORTRAN_FALSE;
}


logical rtdb_open_(const char *filename, const char *mode, integer *handle,
		   const integer flen, const integer mlen)
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

logical rtdb_close_(const integer *handle, const char *mode, const int mlen)
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

logical rtdb_get_info_(const integer *handle, const char *name, 
		       integer *ma_type, integer *nelem, char *date,
		       const int nlen, const int dlen)
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

logical rtdb_put_(const integer *handle, const char *name, const integer *ma_type,
		  const integer *nelem, const void *array, const int nlen)
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

logical rtdb_get_(const integer *handle, const char *name, 
		  const integer *ma_type, const integer *nelem, 
		  void *array, const int nlen)
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

logical rtdb_ma_get_(const integer *handle, const char *name, integer *ma_type,
		     integer *nelem, integer *ma_handle, const int nlen)
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

  if (rtdb_ma_get(hbuf, nbuf, &typebuf, &handbuf, &nelbuf)) {
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

logical rtdb_cput_(const integer *handle, const char *name,
		   const integer *nelem,
		   const char *array, const int nlen, const int alen)
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

  for (i=0, left=sizeof(abuf), next=abuf;
       i<*nelem;
       i++, array+=alen) {
    
    if (!fortchar_to_string(array, alen, next, left)) {
      (void) fprintf(stderr, "rtdb_cput: abuf is too small, need=%d\n", 
		     (int) (alen + sizeof(abuf) - left));
      return FORTRAN_FALSE;
    }

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

logical rtdb_cget_(const integer *handle, const char *name,
		   const integer *nelem,
		   char *array, const int nlen, const int alen)
{
  int hbuf = (int) *handle;
  char nbuf[256];
  char abuf[10240];
  int nelbuf;
  int typebuf;
  int i;
  char *next;

  if (!fortchar_to_string(name, nlen, nbuf, sizeof(nbuf))) {
    (void) fprintf(stderr, "rtdb_cget: nbuf is too small, need=%d\n", 
		   nlen);
    return FORTRAN_FALSE;
  }

  nelbuf = sizeof(abuf);
  typebuf= (int) MT_CHAR;

  if (!rtdb_get(hbuf, nbuf, typebuf, nelbuf, abuf))
    return FORTRAN_FALSE;

  for (i=0, next=strtok(abuf, "\n");
       next;
       i++, array+=alen, next=strtok((char *) 0, "\n")) {

    if (i == *nelem) {
      (void) fprintf(stderr, "rtdb_cget: array has too few elements\n");
      return FORTRAN_FALSE;
    }

    if (!string_to_fortchar(array, alen, next)) {
      (void) fprintf(stderr, "rtdb_cget: array element is too small\n");
      return FORTRAN_FALSE;
    }
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

logical rtdb_delete_(const integer *handle, const char *name, const int nlen)
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

