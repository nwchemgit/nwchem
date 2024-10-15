#include <stdio.h>
#include <string.h>
#include "rtdb.h"
#include "macdecls.h"
#define MXLGTH 32768

#define FORTRAN_TRUE  ((Logical) 1)
#define FORTRAN_FALSE ((Logical) 0)


int fortchar_to_string(const char *f, int flen, char *buf, const int buflen)
{
  while (flen-- && f[flen] == ' ');

  if (flen < 0) flen=0;		/* Empty strings break use of strtok 
				   since consecutive separators are
				   treated as one */

  if ((flen+1) >= buflen) return 0; /* Won't fit */

  flen++;
  buf[flen] = 0;
  while(flen--) buf[flen] = f[flen];

  return 1;
}

int string_to_fortchar(char *f, int flen, char *buf)
{
  int len = (int) strlen(buf), i;

  if (len > flen) return 0; /* Won't fit */

  for (i=0; i<len; i++) f[i] = buf[i];
  for (i=len; i<flen; i++) f[i] = ' ';

  return 1;
}


Logical rtdb_parallel_(const Logical *mode)
{
  /* This causes problems on machines where true != 1 (i.e. intel)
   * so it is better just to pass what we are given
   *
   * int new = (*mode == FORTRAN_TRUE); 
   */
  int new = *mode;
  int old = rtdb_parallel(new);

  if (old) {
    return FORTRAN_TRUE;
  } else {
    return FORTRAN_FALSE;
  }
}

Logical rtdb_open_(const char *filename, const char *mode, Integer *handle,
		   const Integer flen, const Integer mlen)
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
    *handle = (Integer) hbuf;
    return FORTRAN_TRUE;
  } else {
    return FORTRAN_FALSE;
  }
}

Logical rtdb_clone_(const Integer *handle, const char *suffix, const int mlen)
{
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

Logical rtdb_getfname_(const Integer *handle, char *fname, const int mlen)
{
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

Logical rtdb_close_(const Integer *handle, const char *mode, const int mlen)
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

Logical rtdb_get_info_(const Integer *handle, const char *name, 
		       Integer *ma_type, Integer *nelem, char *date,
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

Logical rtdb_put_(const Integer *handle, const char *name, const Integer *ma_type,
		  const Integer *nelem, const void *array, const int nlen)
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

Logical rtdb_get_(const Integer *handle, const char *name, 
		  const Integer *ma_type, const Integer *nelem, 
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

Logical rtdb_ma_get_(const Integer *handle, const char *name, Integer *ma_type,
		     Integer *nelem, Integer *ma_handle, const int nlen)
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
    *ma_type   = (Integer) typebuf;
    *ma_handle = (Integer) handbuf;
    *nelem     = (Integer) nelbuf;

    return FORTRAN_TRUE;
  } else {
    return FORTRAN_FALSE;
  }
}

Logical rtdb_print_(const Integer *handle, const Logical *print_values)
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

Logical rtdb_cput_(const Integer *handle, const char *name,
		   const Integer *nelem,
		   const char *array, const int nlen, const int alen)
{
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
      const char *element = array;
    
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

Logical rtdb_cget_size_(const Integer *handle, const char *name,
		   Integer *nelem,
		   const int nlen, const int alen)
{

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

Logical rtdb_cget_(const Integer *handle, const char *name,
		   const Integer *nelem,
		   char *array, const int nlen, const int alen)
{
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

  if (!rtdb_get(hbuf, nbuf, typebuf, nelbuf, abuf)) {
      return FORTRAN_FALSE;	/* Not there */
  }

  for (i=0, next=strtok(abuf, "\n");
       next;
       i++, array+=alen, next=strtok((char *) 0, "\n")) {
    char *element = array;

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

Logical rtdb_first_(const Integer *handle, char *name, int nlen)
{
  char nbuf[256];
  
  if (rtdb_first((int) *handle, (int) sizeof(nbuf), nbuf) &&
      string_to_fortchar(name, nlen, nbuf))
    return FORTRAN_TRUE;
  else
    return FORTRAN_FALSE;
}

Logical rtdb_next_(const Integer *handle, char *name, int nlen)
{
  char nbuf[256];

  if (rtdb_next((int) *handle, (int) sizeof(nbuf), nbuf) &&
      string_to_fortchar(name, nlen, nbuf)) {
    return FORTRAN_TRUE;
  } else {
    return FORTRAN_FALSE;
  }
}

Logical rtdb_delete_(const Integer *handle, const char *name, const int nlen)
{
  int hbuf = (int) *handle;
  char nbuf[256];

  if (!fortchar_to_string(name, nlen, nbuf, sizeof(nbuf))) {
    (void) fprintf(stderr, "rtdb_delete: nbuf is too small, need=%d\n", 
		   nlen);
    return FORTRAN_FALSE;
  }

  if (rtdb_delete(hbuf, nbuf)) {
    return FORTRAN_TRUE;
  } else {
    return FORTRAN_FALSE;
  }
}

extern void rtdb_print_usage(void);

void rtdb_print_usage_()
{
  rtdb_print_usage();
}
