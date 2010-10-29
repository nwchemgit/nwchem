/*$Id$*/
#include <stdio.h>
#include <string.h>
#include "srtdb.h"
#include "macdecls.h"
#ifdef CRAY
#include <fortran.h>
#endif

#define FORTRAN_TRUE  ((Logical) 1)
#define FORTRAN_FALSE ((Logical) 0)


#if (defined(CRAY) || defined(USE_FCD)) && !defined(__crayx1)
Logical FATR srtdb_open_(_fcd filename, _fcd mode, Integer *handle)
{
  int flen = _fcdlen(filename);
  int mlen = _fcdlen(mode);
#else
Logical FATR srtdb_open_(const char *filename, const char *mode, Integer *handle,
		   const Integer flen, const Integer mlen)
{
#endif
  char fbuf[256], mbuf[256];
  int hbuf;

  if (!fortchar_to_string(filename, flen, fbuf, sizeof(fbuf))) {
    (void) fprintf(stderr, "srtdb_open: fbuf is too small, need=%d\n",
		   (int) flen);
    return FORTRAN_FALSE;
  }

  if (!fortchar_to_string(mode, mlen, mbuf, sizeof(mbuf))) {
    (void) fprintf(stderr, "srtdb_open: mbuf is too small, need=%d\n",
		   (int) mlen);
    return FORTRAN_FALSE;
  }

  if (srtdb_open(fbuf, mbuf, &hbuf)) {
    *handle = (Integer) hbuf;
    return FORTRAN_TRUE;
  }
  else {
    return FORTRAN_FALSE;
  }
}


#if (defined(CRAY) || defined(USE_FCD)) && !defined(__crayx1)
Logical FATR srtdb_close_(const Integer *handle, _fcd mode)
{
  int mlen = _fcdlen(mode);
#else
Logical FATR srtdb_close_(const Integer *handle, const char *mode, const int mlen)
{
#endif
  char mbuf[256];
  int hbuf = (int) *handle;

  if (!fortchar_to_string(mode, mlen, mbuf, sizeof(mbuf))) {
    (void) fprintf(stderr, "srtdb_close: mbuf is too small, need=%d\n", mlen);
    return FORTRAN_FALSE;
  }
 if (srtdb_close(hbuf, mbuf))
    return FORTRAN_TRUE;
 else
    return FORTRAN_FALSE;
}
#if (defined(CRAY) || defined(USE_FCD)) && !defined(__crayx1)
Logical FATR srtdb_get_info_(const Integer *handle, _fcd name, 
		       Integer *ma_type, Integer *nelem, _fcd date)
{
    int nlen = _fcdlen(name);
    int dlen = _fcdlen(date);
#else
Logical FATR srtdb_get_info_(const Integer *handle, const char *name, 
		       Integer *ma_type, Integer *nelem, char *date,
		       const int nlen, const int dlen)
{
#endif

  int hbuf = (int) *handle;
  char dbuf[26], nbuf[256];
  int nelbuf, typebuf;

  if (!fortchar_to_string(name, nlen, nbuf, sizeof(nbuf))) {
    (void) fprintf(stderr, "srtdb_get_info: nbuf is too small, need=%d\n", 
		   nlen);
    return FORTRAN_FALSE;
  }

  if (dlen < 24) {
    (void) fprintf(stderr, "srtdb_get_info: date must be > character*24\n");
    return FORTRAN_FALSE;
  }
    
  if (srtdb_get_info(hbuf, nbuf, &typebuf, &nelbuf, dbuf)) {
    *ma_type = (Integer) typebuf;
    *nelem   = (Integer) nelbuf;

    if (typebuf == MT_CHAR)	/* Fortran is ignorant of trailing null char */
      *nelem = *nelem - 1;

    if (!string_to_fortchar(date, dlen, dbuf)) {
      (void) fprintf(stderr, "srtdb_get_info: nbuf is too small, need=%d\n", 
		     nlen);
      return FORTRAN_FALSE;
    }

    return FORTRAN_TRUE;
  }
  else {
    return FORTRAN_FALSE;
  }
}
#if (defined(CRAY) || defined(USE_FCD)) && !defined(__crayx1)
Logical FATR srtdb_put_(const Integer *handle, _fcd name, const Integer *ma_type,
		  const Integer *nelem, const void *array)
{
    int nlen = _fcdlen(name);
#else
Logical FATR srtdb_put_(const Integer *handle, const char *name, const Integer *ma_type,
		  const Integer *nelem, const void *array, const int nlen)
{
#endif
  int hbuf = (int) *handle;
  char nbuf[256];
  int nelbuf;
  int typebuf;

  if (!fortchar_to_string(name, nlen, nbuf, sizeof(nbuf))) {
    (void) fprintf(stderr, "srtdb_put: nbuf is too small, need=%d\n", 
		   nlen);
    return FORTRAN_FALSE;
  }

  nelbuf = (int) *nelem;
  typebuf= (int) *ma_type;

#ifdef DEBUG
  printf("put: rtdb=%d, mat=%d, nel=%d, name=%s\n", hbuf, typebuf, nelbuf, nbuf);
  fflush(stdout);
#endif

  if (srtdb_put(hbuf, nbuf, typebuf, nelbuf, array))
    return FORTRAN_TRUE;
  else
    return FORTRAN_FALSE;
}
#if (defined(CRAY) || defined(USE_FCD)) && !defined(__crayx1)
Logical FATR srtdb_get_(const Integer *handle, _fcd name, 
		  const Integer *ma_type, const Integer *nelem, 
		  void *array)
{
    int nlen = _fcdlen(name);
#else
Logical FATR srtdb_get_(const Integer *handle, const char *name, 
		  const Integer *ma_type, const Integer *nelem, 
		  void *array, const int nlen)
{
#endif
  int hbuf = (int) *handle;
  char nbuf[256];
  int nelbuf;
  int typebuf;

  if (!fortchar_to_string(name, nlen, nbuf, sizeof(nbuf))) {
    (void) fprintf(stderr, "srtdb_get: nbuf is too small, need=%d\n", 
		   nlen);
    return FORTRAN_FALSE;
  }

  nelbuf = (int) *nelem;
  typebuf= (int) *ma_type;

#ifdef DEBUG
  printf("get: rtdb=%d, mat=%d, nel=%d, name=%s\n", hbuf, typebuf, nelbuf, nbuf);
  fflush(stdout);
#endif

  if (srtdb_get(hbuf, nbuf, typebuf, nelbuf, array)) {
    return FORTRAN_TRUE;
  }
  else
    return FORTRAN_FALSE;
}


Logical FATR srtdb_print_(const Integer *handle, const Logical *print_values)
{
  int hbuf = (int) *handle;
  int pbuf = (int) *print_values;

  if (srtdb_print(hbuf, pbuf))
    return FORTRAN_TRUE;
  else
    return FORTRAN_FALSE;
}

#if (defined(CRAY) || defined(USE_FCD)) && !defined(__crayx1)
Logical FATR srtdb_cput_(const Integer *handle, _fcd name,
		   const Integer *nelem,
		   _fcd farray)
{
    int nlen = _fcdlen(name);
    int alen = _fcdlen(farray);
    char *array = _fcdtocp(farray);
#else
Logical FATR srtdb_cput_(const Integer *handle, const char *name,
		   const Integer *nelem,
		   const char *array, const int nlen, const int alen)
{
#endif
/*
  Insert an array of Fortran character variables into the data base.
  Each array element is striped of trailing blanks, terminated with CR,
  and appended to the list. The entire array must fit into abuf.
*/

  int hbuf = (int) *handle;
  char nbuf[256];
  char abuf[20480]=" ";
  int nelbuf;
  int typebuf;
  int i, left;
  char *next;
  for (i=0, left=sizeof(abuf), next=abuf;
       i<*nelem;
       i++, array+=alen) {
#if defined(CRAY) && !defined(__crayx1)
      _fcd element = _cptofcd(array, alen);
#elif defined(WIN32)
      _fcd element;
      element.string = array;
      element.len = alen;
#elif defined(USE_FCD)
#error Do something about _fcd
#else
      const char *element = array;
#endif
    
    if (!fortchar_to_string(element, alen, next, left)) {
      (void) fprintf(stderr, "srtdb_cput: abuf is too small, need=%d\n", 
		     (int) (alen + sizeof(abuf) - left));
      return FORTRAN_FALSE;
    }
    left -= strlen(next) + 1;
    next += strlen(next) + 1;
    if (i != (*nelem - 1))
      *(next-1) = '\n';
  }

  if (!fortchar_to_string(name, nlen, nbuf, sizeof(nbuf))) {
    (void) fprintf(stderr, "srtdb_cput: nbuf is too small, need=%d\n", 
		   nlen);
    return FORTRAN_FALSE;
  }

  nelbuf = strlen(abuf) + 1;
  typebuf= (int) MT_CHAR;

#ifdef DEBUG
  printf("cput: rtdb=%d, mat=%d, nel=%d, name=%s\n", hbuf, typebuf, nelbuf, nbuf);
  fflush(stdout);
#endif

  if (srtdb_put(hbuf, nbuf, typebuf, nelbuf, abuf))
    return FORTRAN_TRUE;
  else
    return FORTRAN_FALSE;
}


#if (defined(CRAY) || defined(USE_FCD)) && !defined(__crayx1)
Logical FATR srtdb_cget_(const Integer *handle, _fcd name,
		   const Integer *nelem,
		   _fcd farray)
{
    int nlen = _fcdlen(name);
    int alen = _fcdlen(farray);
    char *array = _fcdtocp(farray);
#else
Logical FATR srtdb_cget_(const Integer *handle, const char *name,
		   const Integer *nelem,
		   char *array, const int nlen, const int alen)
{
#endif
/*
  Read an array of Fortran character variables from the data base.

  Put stored the array as follows:
  .  Each array element is striped of trailing blanks, terminated with CR,
  .  and appended to the list. The entire array must fit into abuf.
*/

  int hbuf = (int) *handle;
  char nbuf[256];
  char abuf[20480];
  /*  char abuf[10240];*/
  int nelbuf;
  int typebuf;
  int i;
  char *next;

  if (!fortchar_to_string(name, nlen, nbuf, sizeof(nbuf))) {
    (void) fprintf(stderr, "srtdb_cget: nbuf is too small, need=%d\n", 
		   nlen);
    return FORTRAN_FALSE;
  }

  nelbuf = sizeof(abuf);
  typebuf= (int) MT_CHAR;

#ifdef DEBUG
  printf("cget: rtdb=%d, mat=%d, nel=%d, name=%s\n", hbuf, typebuf, nelbuf, nbuf);
  fflush(stdout);
#endif

  if (!srtdb_get(hbuf, nbuf, typebuf, nelbuf, abuf))
      return FORTRAN_FALSE;	/* Not there */

  for (i=0, next=strtok(abuf, "\n");
       next;
       i++, array+=alen, next=strtok((char *) 0, "\n")) {
#if defined(CRAY) && !defined(__crayx1)
      _fcd element = _cptofcd(array, alen);
#elif defined(WIN32)
      _fcd element;
      element.string = array;
      element.len = alen;
#elif defined(USE_FCD)
#error Do something about _fcd
#else
      char *element = array;
#endif

    if (i == *nelem) {
      (void) fprintf(stderr, "srtdb_cget: array has too few elements\n");
      (void) fprintf(stderr, "srtdb_cget: name was <<%s>>\n",name);
      return FORTRAN_FALSE;
    }

    if (!string_to_fortchar(element, alen, next)) {
      (void) fprintf(stderr, "srtdb_cget: array element is too small\n");
      (void) fprintf(stderr, "srtdb_cget: name was <<%s>>\n",name);
      return FORTRAN_FALSE;
    }
  }
  return FORTRAN_TRUE;
}


#if (defined(_CRAY) || defined(USE_FCD)) && !defined(__crayx1)
Logical FATR srtdb_first_(const Integer *handle, _fcd name)
#else
Logical FATR srtdb_first_(const Integer *handle, char *name, int nlen)
#endif
{
#if (defined(_CRAY) || defined(USE_FCD)) && !defined(__crayx1)
  // dummy arg, value reassigned by string_to_fortchar in this case
  int nlen = _fcdlen(name);
#endif
  char nbuf[256];
  
  if (srtdb_first((int) *handle, (int) sizeof(nbuf), nbuf) &&
      string_to_fortchar(name, nlen, nbuf))
    return FORTRAN_TRUE;
  else
    return FORTRAN_FALSE;
}

#if (defined(_CRAY) || defined(USE_FCD)) && !defined(__crayx1)
Logical FATR srtdb_next_(const Integer *handle, _fcd name)
#else
Logical FATR srtdb_next_(const Integer *handle, char *name, int nlen)
#endif
{
#if (defined(_CRAY) || defined(USE_FCD)) && !defined(__crayx1)
  // dummy arg, value reassigned by string_to_fortchar in this case
  int nlen = _fcdlen(name);
#endif
  char nbuf[256];

  if (srtdb_next((int) *handle, (int) sizeof(nbuf), nbuf) &&
      string_to_fortchar(name, nlen, nbuf))
    return FORTRAN_TRUE;
  else
    return FORTRAN_FALSE;
}
#if (defined(CRAY) || defined(USE_FCD)) && !defined(__crayx1)
Logical FATR srtdb_delete_(const Integer *handle, _fcd name)
{
  int nlen = _fcdlen(name);
#else
Logical FATR srtdb_delete_(const Integer *handle, const char *name, const int nlen)
{
#endif
  int hbuf = (int) *handle;
  char nbuf[256];

  if (!fortchar_to_string(name, nlen, nbuf, sizeof(nbuf))) {
    (void) fprintf(stderr, "srtdb_delete: nbuf is too small, need=%d\n", 
		   nlen);
    return FORTRAN_FALSE;
  }

  if (srtdb_delete(hbuf, nbuf))
    return FORTRAN_TRUE;
  else
    return FORTRAN_FALSE;
}

extern void srtdb_print_usage(void);

void FATR srtdb_print_usage_()
{
  srtdb_print_usage();
}
