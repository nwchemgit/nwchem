/*$Id$*/
#include <stdio.h>
#include <string.h>
#include "rtdb.h"
#include "macdecls.h"
#ifdef CRAY
#include <fortran.h>
#endif
#define FATR 
#define MXLGTH 32768

#define FORTRAN_TRUE  ((Logical) 1)
#define FORTRAN_FALSE ((Logical) 0)


#if (defined(CRAY) || defined(USE_FCD)) && !defined(__crayx1)
int fortchar_to_string(_fcd f, int flen, char *buf, const int buflen)
#else
int fortchar_to_string(const char *f, int flen, char *buf, const int buflen)
#endif
{
#if (defined(CRAY) || defined(USE_FCD)) && !defined(__crayx1)
  char *fstring;
    fstring = _fcdtocp(f);
    flen = _fcdlen(f);

  while (flen-- && fstring[flen] == ' ')
    ;

  if (flen < 0) flen=0;		/* Empty strings break use of strtok 
				   since consecutive separators are
				   treated as one */
  if ((flen+1) >= buflen)
    return 0;			/* Won't fit */

  flen++;
  buf[flen] = 0;
  while(flen--)
    buf[flen] = fstring[flen];
#else 
  while (flen-- && f[flen] == ' ')
    ;

  if (flen < 0) flen=0;		/* Empty strings break use of strtok 
				   since consecutive separators are
				   treated as one */
  if ((flen+1) >= buflen)
    return 0;			/* Won't fit */

  flen++;
  buf[flen] = 0;
  while(flen--)
    buf[flen] = f[flen];
#endif 

  return 1;
}

#if (defined(CRAY) || defined(USE_FCD)) && !defined(__crayx1)
int string_to_fortchar(_fcd f, int flen, char *buf)
#else
int string_to_fortchar(char *f, int flen, char *buf)
#endif
{
  int len = (int) strlen(buf), i;
#if (defined(CRAY) || defined(USE_FCD)) && !defined(__crayx1)
  flen = _fcdlen(f);
#endif

  if (len > flen) {
    return 0;			/* Won't fit */
}
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


Logical FATR rtdb_parallel_(const Logical *mode)
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

#if (defined(CRAY) || defined(USE_FCD)) && !defined(__crayx1)
Logical FATR rtdb_open_(_fcd filename, _fcd mode, Integer *handle)
{
  int flen = _fcdlen(filename);
  int mlen = _fcdlen(mode);
#else
Logical FATR rtdb_open_(const char *filename, const char *mode, Integer *handle,
		   const Integer flen, const Integer mlen)
{
#endif
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
    *handle = (Integer) hbuf;
    return FORTRAN_TRUE;
  }
  else {
    return FORTRAN_FALSE;
  }
}
#if (defined(CRAY) || defined(USE_FCD)) && !defined(__crayx1)
Logical FATR rtdb_clone_(const Integer *handle, _fcd suffix)
{
  int mlen = _fcdlen(suffix);
#else
Logical FATR rtdb_clone_(const Integer *handle, const char *suffix, const int mlen)
{
#endif
  char mbuf[256];
  int hbuf = (int) *handle;

  if (!fortchar_to_string(suffix, mlen, mbuf, sizeof(mbuf))) {
    (void) fprintf(stderr, "rtdb_clone: mbuf is too small, need=%d\n", mlen);
    return FORTRAN_FALSE;
  }
 if (rtdb_clone(hbuf, mbuf))
    return FORTRAN_TRUE;
 else
    return FORTRAN_FALSE;
}
#if (defined(CRAY) || defined(USE_FCD)) && !defined(__crayx1)
Logical FATR rtdb_getfname_(const Integer *handle, _fcd fname)
{
  int mlen = _fcdlen(fname);
#else
Logical FATR rtdb_getfname_(const Integer *handle, char *fname, const int mlen)
{
#endif
  char mbuf[256];
  int hbuf = (int) *handle;

  if (rtdb_getfname(hbuf, mbuf)){
    if (!string_to_fortchar(fname, mlen, mbuf)) {
      (void) fprintf(stderr, "rtdb_fname: mbuf is too small, need=%ld\n",strlen(mbuf));
      return FORTRAN_FALSE;
    }
    return FORTRAN_TRUE;
  }
 else
    return FORTRAN_FALSE;
}
#if (defined(CRAY) || defined(USE_FCD)) && !defined(__crayx1)
Logical FATR rtdb_close_(const Integer *handle, _fcd mode)
{
  int mlen = _fcdlen(mode);
#else
Logical FATR rtdb_close_(const Integer *handle, const char *mode, const int mlen)
{
#endif
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
#if (defined(CRAY) || defined(USE_FCD)) && !defined(__crayx1)
Logical FATR rtdb_get_info_(const Integer *handle, _fcd name, 
		       Integer *ma_type, Integer *nelem, _fcd date)
{
    int nlen = _fcdlen(name);
    int dlen = _fcdlen(date);
#else
Logical FATR rtdb_get_info_(const Integer *handle, const char *name, 
		       Integer *ma_type, Integer *nelem, char *date,
		       const int nlen, const int dlen)
{
#endif

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
    *ma_type = (Integer) typebuf;
    *nelem   = (Integer) nelbuf;

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
#if (defined(CRAY) || defined(USE_FCD)) && !defined(__crayx1)
Logical FATR rtdb_put_(const Integer *handle, _fcd name, const Integer *ma_type,
		  const Integer *nelem, const void *array)
{
    int nlen = _fcdlen(name);
#else
Logical FATR rtdb_put_(const Integer *handle, const char *name, const Integer *ma_type,
		  const Integer *nelem, const void *array, const int nlen)
{
#endif
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
#if (defined(CRAY) || defined(USE_FCD)) && !defined(__crayx1)
Logical FATR rtdb_get_(const Integer *handle, _fcd name, 
		  const Integer *ma_type, const Integer *nelem, 
		  void *array)
{
    int nlen = _fcdlen(name);
#else
Logical FATR rtdb_get_(const Integer *handle, const char *name, 
		  const Integer *ma_type, const Integer *nelem, 
		  void *array, const int nlen)
{
#endif
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
#if (defined(CRAY) || defined(USE_FCD)) && !defined(__crayx1)
Logical FATR rtdb_ma_get_(const Integer *handle, _fcd name, Integer *ma_type,
		     Integer *nelem, Integer *ma_handle)
{
    int nlen = _fcdlen(name);
#else
Logical FATR rtdb_ma_get_(const Integer *handle, const char *name, Integer *ma_type,
		     Integer *nelem, Integer *ma_handle, const int nlen)
{
#endif
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
    *ma_type   = (Integer) typebuf;
    *ma_handle = (Integer) handbuf;
    *nelem     = (Integer) nelbuf;

    return FORTRAN_TRUE;
  }
  else
    return FORTRAN_FALSE;
}

Logical FATR rtdb_print_(const Integer *handle, const Logical *print_values)
{
  int hbuf = (int) *handle;
  int pbuf = (int) *print_values;

  if (rtdb_print(hbuf, pbuf))
    return FORTRAN_TRUE;
  else
    return FORTRAN_FALSE;
}

/**
\ingroup rtdb
@{
*/

/** 
  \brief Store a character string on the RTDB

  This function is supposed to be called from a Fortran code.

  \param handle [Input] the RTDB handle
  \param name   [Input] the key
  \param nelem  [Input] the length of the character buffer
  \param array  [Input] the value of the string

  \return Return FORTRAN_TRUE if successfull, and FORTRAN_FALSE otherwise.
*/
#if (defined(CRAY) || defined(USE_FCD)) && !defined(__crayx1)
Logical FATR rtdb_cput_(const Integer *handle, _fcd name,
		   const Integer *nelem,
		   _fcd farray)
{
    int nlen = _fcdlen(name);
    int alen = _fcdlen(farray);
    char *array = _fcdtocp(farray);
#else
Logical FATR rtdb_cput_(const Integer *handle, const char *name,
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
  char abuf[MXLGTH]=" ";
  int nelbuf; 
  int typebuf;
  int i, left;
  char *next;

  for (i=0, left=sizeof(abuf), next=abuf;
       i<*nelem;
       i++, array+=alen) {
#if defined(CRAY) && !defined(__crayx1)
      _fcd element = _cptofcd(array, alen);
#elif defined(WIN32) &&! defined(__MINGW32__)
      _fcd element;
      element.string = array;
      element.len = alen;
#elif defined(USE_FCD)
#error Do something about _fcd
#else
      const char *element = array;
#endif
    
    if (!fortchar_to_string(element, alen, next, left)) {
      (void) fprintf(stderr, "rtdb_cput: abuf is too small, increase MXLGTH to=%d\n", 
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

/*MV*/
/** 
  \brief Retrieve the length of a character string on the RTDB

  When a character string is stored on the RTDB it is not necessarily obvious
  how long the string is. Hence it is useful to be able to check the length of
  the string before attempting to retrieve it to make sure that a buffer of
  sufficient length is provided.

  This function is supposed to be called from a Fortran code.

  \param handle [Input] the RTDB handle
  \param name   [Input] the key
  \param nelem  [Output] the number of characters

  \return Return FORTRAN_TRUE if successfull, and FORTRAN_FALSE otherwise.
*/
#if (defined(CRAY) || defined(USE_FCD)) && !defined(__crayx1)
Logical FATR rtdb_cget_size_(const Integer *handle, _fcd name,
		   const Integer *nelem)
{
    int nlen = _fcdlen(name);
    int alen = _fcdlen(farray);
#else
Logical FATR rtdb_cget_size_(const Integer *handle, const char *name,
		   Integer *nelem,
		   const int nlen, const int alen)
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
  char abuf[MXLGTH];
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

  if (!rtdb_get(hbuf, nbuf, typebuf, nelbuf, abuf)) {
      return FORTRAN_FALSE;	/* Not there */
  }

  for (i=0, next=strtok(abuf, "\n");
       next;
       i++, next=strtok((char *) 0, "\n")) {
  }
  *nelem = i;
  return FORTRAN_TRUE;
}




/*MV*/
/** 
  \brief Retrieve a character string from the RTDB

  Retrieve a character string from the RTDB and return it in the character
  buffer provided. The buffer must be large enough to hold the value otherwise
  the function will fail.

  This function is supposed to be called from a Fortran code.

  \param handle [Input] the RTDB handle
  \param name   [Input] the key
  \param nelem  [Input] the length of the character buffer
  \param array  [Output] the value of the string

  \return Return FORTRAN_TRUE if successfull, and FORTRAN_FALSE otherwise.
*/
#if (defined(CRAY) || defined(USE_FCD)) && !defined(__crayx1)
Logical FATR rtdb_cget_(const Integer *handle, _fcd name,
		   const Integer *nelem,
		   _fcd farray)
{
    int nlen = _fcdlen(name);
    int alen = _fcdlen(farray);
    char *array = _fcdtocp(farray);
#else
Logical FATR rtdb_cget_(const Integer *handle, const char *name,
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
  char abuf[MXLGTH];
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

#ifdef DEBUG
  printf("cget: rtdb=%d, mat=%d, nel=%d, name=%s\n", hbuf, typebuf, nelbuf, nbuf);
  fflush(stdout);
#endif

  if (!rtdb_get(hbuf, nbuf, typebuf, nelbuf, abuf))
      return FORTRAN_FALSE;	/* Not there */

  for (i=0, next=strtok(abuf, "\n");
       next;
       i++, array+=alen, next=strtok((char *) 0, "\n")) {
#if defined(CRAY) && !defined(__crayx1)
      _fcd element = _cptofcd(array, alen);
#elif defined(WIN32) &&! defined(__MINGW32__)
      _fcd element;
      element.string = array;
      element.len = alen;
#elif defined(USE_FCD)
#error Do something about _fcd
#else
      char *element = array;
#endif

    if (i == *nelem) {
      (void) fprintf(stderr, "rtdb_cget: array has too few elements\n");
      (void) fprintf(stderr, "rtdb_cget: name was <<%s>>\n",name);
      return FORTRAN_FALSE;
    }

    if (!string_to_fortchar(element, alen, next)) {
      (void) fprintf(stderr, "rtdb_cget: array element is too small\n");
      (void) fprintf(stderr, "rtdb_cget: name was <<%s>>\n",name);
      return FORTRAN_FALSE;
    }
  }
  return FORTRAN_TRUE;
}

/**
@}
*/

#if (defined(_CRAY) || defined(USE_FCD)) && !defined(__crayx1)
Logical FATR rtdb_first_(const Integer *handle, _fcd name)
#else
Logical FATR rtdb_first_(const Integer *handle, char *name, int nlen)
#endif
{
#if (defined(_CRAY) || defined(USE_FCD)) && !defined(__crayx1)
  // dummy arg, value reassigned by string_to_fortchar in this case
  int nlen = _fcdlen(name);
#endif
  char nbuf[256];
  
  if (rtdb_first((int) *handle, (int) sizeof(nbuf), nbuf) &&
      string_to_fortchar(name, nlen, nbuf))
    return FORTRAN_TRUE;
  else
    return FORTRAN_FALSE;
}

#if (defined(_CRAY) || defined(USE_FCD)) && !defined(__crayx1)
Logical FATR rtdb_next_(const Integer *handle, _fcd name)
#else
Logical FATR rtdb_next_(const Integer *handle, char *name, int nlen)
#endif
{
#if (defined(_CRAY) || defined(USE_FCD)) && !defined(__crayx1)
  // dummy arg, value reassigned by string_to_fortchar in this case
  int nlen = _fcdlen(name);
#endif
  char nbuf[256];

  if (rtdb_next((int) *handle, (int) sizeof(nbuf), nbuf) &&
      string_to_fortchar(name, nlen, nbuf))
    return FORTRAN_TRUE;
  else
    return FORTRAN_FALSE;
}
#if (defined(CRAY) || defined(USE_FCD)) && !defined(__crayx1)
Logical FATR rtdb_delete_(const Integer *handle, _fcd name)
{
  int nlen = _fcdlen(name);
#else
Logical FATR rtdb_delete_(const Integer *handle, const char *name, const int nlen)
{
#endif
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

extern void rtdb_print_usage(void);

void FATR rtdb_print_usage_()
{
  rtdb_print_usage();
}
