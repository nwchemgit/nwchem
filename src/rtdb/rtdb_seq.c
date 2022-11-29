/*$Id$*/
#if defined(MACX)
#include <stdio.h>
#endif
#include <stdlib.h>
#include <sys/types.h>
#if defined(CATAMOUNT)
#warning including iobufffff
#include </opt/xt-libc/default/amd64/include/stdio.h>
#define USE_IOBUF_MACROS
#endif
#include <stdio.h>

#ifdef USE_HDBM
#include "hdbm/hdbm.h"
#define DBT datum
#define SIZE(a) ((a).dsize)
#define DATA(a) ((a).dptr)
#define WRAP(data, size, p) datum_wrap(data, size, p)
#endif
#ifdef USE_DB
#include <db.h>
#define SIZE(a) ((a).size)
#define DATA(a) ((a).data)
#define WRAP(data, size, p) *(p) = wrap_DBT(data, size)
#endif

#ifdef WIN32
#include <io.h>
#define R_OK 4
#define W_OK 2
#else
#include <unistd.h>
#endif
#include <string.h>
#include <time.h>
#if !defined(IPSC) && !defined(WIN32)
#include <sys/time.h>
#endif
#include <sys/stat.h>
#include <fcntl.h>
#include "rtdb.h"
#include "macdecls.h"
#include "misc.h"

#if !defined(LINUX)
extern char *strdup(const char *);
extern void *malloc(size_t);
extern void free(void *);
#endif

#define MAX_RTDB 100

#define FORTRAN_TRUE  1L
#define FORTRAN_FALSE 0L

static struct {			/* Keep track of active RTDBs */
  int active;
  char *filename;
  int scratch;
#ifdef USE_HDBM
  hdbm db;
#endif
#ifdef USE_DB
  DB *db;
#endif
} rtdb[MAX_RTDB];

struct info_struct{		/* Matching info entry for each datum */
  int ma_type;
  int nelem;
  char date[26];
};

static char info_header[] = "!rtdb!"; /* Prefix for info entries */

#ifdef USE_DB
static DBT wrap_DBT(const void *data, size_t size)
/*
  Make a DBT that refers to the region pointed to
  by the arguments .. NO MEMORY IS ALLOCATED
*/
{
  DBT tmp;

  tmp.data = (void *) data; tmp.size = size;

  return tmp;
}
#endif

static long file_size(const char *filename)
/*
  Return size of file in bytes or -1 if does not exist
*/
{
  FILE *file = fopen(filename, "r+");
  
  long length;

  if (!file)
    return (size_t) -1;

  (void) fseek(file, 0L, 2);
  length = ftell(file);
  (void) fclose(file);

  /* printf("file %s is %ld bytes long\n", filename, length); */

  return length;
}

void ma_print(FILE *file, const int ma_type, const int nelem, void *p)
{
  int i, nprint;

  switch (ma_type) {
  case MT_C_CHAR:	/* char */
  case MT_F_BYTE:	/* Fortran byte */

    (void) fprintf(file, "%.*s\n", nelem, (char *) p);
    break;

  case MT_C_INT:	/* int */
    for (nprint=i=0; i<nelem; i++) {
      nprint += fprintf(file, "%d ", ((int *) p)[i]);
      if (nprint >= 72) {
	(void) fprintf(file, "\n");
	nprint = 0;
      }
    }
    if (nprint > 0) (void) fprintf(file, "\n");
    break;

  case MT_F_INT:	/* Fortran integer ... not equivalent on KSR */
    for (nprint=i=0; i<nelem; i++) {
      nprint += fprintf(file, "%ld ", ((Integer *) p)[i]);
      if (nprint >= 72) {
	(void) fprintf(file, "\n");
	nprint = 0;
      }
    }
    if (nprint > 0) (void) fprintf(file, "\n");
    break;


  case MT_C_LONGINT:	/* long int */

    for (nprint=i=0; i<nelem; i++) {
      nprint += fprintf(file, "%ld ", ((long *) p)[i]);
      if (nprint >= 72) {
	(void) fprintf(file, "\n");
	nprint = 0;
      }
    }
    if (nprint > 0) (void) fprintf(file, "\n");
    break;

  case MT_C_FLOAT:	/* float */
  case MT_F_REAL:	/* Fortran real */

    for (nprint=i=0; i<nelem; i++) {
      nprint += fprintf(file, "%.7e ", ((float *) p)[i]);
      if (nprint >= 72) {
	(void) fprintf(file, "\n");
	nprint = 0;
      }
    }
    if (nprint > 0) (void) fprintf(file, "\n");
    break;

  case MT_C_DBL:	/* double */
  case MT_F_DBL:	/* Fortran double precision */

    for (nprint=i=0; i<nelem; i++) {
      nprint += fprintf(file, "%.14e ", ((double *) p)[i]);
      if (nprint >= 72) {
	(void) fprintf(file, "\n");
	nprint = 0;
      }
    }
    if (nprint > 0) (void) fprintf(file, "\n");
    break;

  case MT_C_SCPL:	/* single precision complex */
  case MT_F_SCPL:	/* Fortran single precision complex */

    for (nprint=i=0; i<nelem; i++) {
      nprint += fprintf(file, "(%.7e,%.7e) ", ((float *) p)[2*i], 
			((float *) p)[2*i+1]);
      if (nprint >= 72) {
	(void) fprintf(file, "\n");
	nprint = 0;
      }
    }
    if (nprint > 0) (void) fprintf(file, "\n");
    break;

  case MT_C_DCPL:	/* double precision complex */
  case MT_F_DCPL:	/* Fortran double precision complex */

    for (nprint=i=0; i<nelem; i++) {
      nprint += fprintf(file, "(%.14e,%.14e) ", ((double *) p)[2*i], 
			((double *) p)[2*i+1]);
      if (nprint >= 72) {
	(void) fprintf(file, "\n");
	nprint = 0;
      }
    }
    if (nprint > 0) (void) fprintf(file, "\n");
    break;

  case MT_F_LOG:	/* Fortran logical */
    for (nprint=i=0; i<nelem; i++) {
      if (((int *) p)[i] == FORTRAN_TRUE)
	nprint += fprintf(file, "t ");
      else
	nprint += fprintf(file, "f ");
      if (nprint >= 72) {
	(void) fprintf(file, "\n");
	nprint = 0;
      }
    }
    if (nprint > 0) (void) fprintf(file, "\n");
    break;


  case MT_C_LDCPL:	/* long double precision complex */
  case MT_C_LDBL:	/* long double */

    (void) fprintf(file, " !! printing long double not supported !!\n");
    break;

  default:

    (void) fprintf(file, " !! %d = unknown data type\n", ma_type);
  }
}

const char *ma_typename(const int ma_type)
{
  switch (ma_type) {
  case MT_C_CHAR:	/* char */
  case MT_F_BYTE:	/* Fortran byte */
    return "char"; break;

  case MT_C_INT:	/* int */
  case MT_F_INT:	/* Fortran integer ... not equivalent on KSR */

    return "int"; break;

  case MT_C_LONGINT:	/* long int */

    return "long"; break;

  case MT_C_FLOAT:	/* float */
  case MT_F_REAL:	/* Fortran real */

    return "float"; break;

  case MT_C_DBL:	/* double */
  case MT_F_DBL:	/* Fortran double precision */

    return "double"; break;

  case MT_C_SCPL:	/* single precision complex */
  case MT_F_SCPL:	/* Fortran single precision complex */

    return "complex float";

  case MT_C_DCPL:	/* double precision complex */
  case MT_F_DCPL:	/* Fortran double precision complex */

    return "complex double"; break;

  case MT_F_LOG:	/* Fortran logical */

    return "logical"; break;

  case MT_C_LDCPL:	/* long double precision complex */
  case MT_C_LDBL:	/* long double */
  default:

    return "invalid"; break;
  }
}

void rtdb_init_()
{
  int new;
  
  for (new=0; new<MAX_RTDB; new++)
    rtdb[new].active = 0;

}

int rtdb_seq_open(const char *filename, const char *mode, int *handle)
/*
    Filename = path to file associated with the data base
    mode     = 'new'     Open only if it does not exist already
               'old',    Open only if it does exist already
               'unknown' Create new or open existing (preserving contents)
               'empty'   Create new or open existing (deleting contents)
               'scratch' Create new or open existing (deleting contents)
                         and automatically delete upon closing.  Also, items
                         cached in memory are not written to disk.

    handle   = returns handle by which all future references to the
               data base are made
*/
{
#ifdef USE_DB
  HASHINFO openinfo;
#endif
  int flags = O_RDWR | O_CREAT;
#if defined(__bgq__)
  int exists = access(filename, F_OK) == 0;
#else
  int exists = access(filename, R_OK | W_OK) == 0;
#endif
  int new;

  /* See if the data base is already open ... if so return and complain */

  for (new=0; new<MAX_RTDB; new++)
    if (rtdb[new].active)
      if (strcmp(filename, rtdb[new].filename) == 0) {
	(void) fprintf(stderr, "rtdb_seq_open: %s is already open\n", 
		       filename);
	return 0;
      }

  /* Figure out if the file exists */

  if (strcmp(mode, "new") == 0) {
    if (exists) {
      (void) fprintf(stderr, 
		     "rtdb_seq_open: %s exists but is being opened new\n",
		     filename);
      return 0;
    }
  }
  else if (strcmp(mode, "old") == 0) {
    if (!exists) {
      (void) fprintf(stderr, 
		     "rtdb_seq_open: %s does not exist, cannot open old\n",
		     filename);
      return 0;
    }
  }
  else if (strcmp(mode, "empty") == 0) {
    flags |= O_TRUNC;
  }
  else if (strcmp(mode, "scratch") == 0) {
    flags |= O_TRUNC;
  }
  else if (strcmp(mode, "unknown") == 0)
    ;
  else {
    (void) fprintf(stderr, "rtdb_seq_open: unknown mode=%s\n", mode);
    return 0;
  }

  /* Find next free rtdb entry */

  for (new=0; new<MAX_RTDB && rtdb[new].active; new++)
    ;
  if (new == MAX_RTDB) {
    (void) fprintf(stderr, "rtdb_seq_open: too many data bases open, max=%d\n",
		   MAX_RTDB);
    return 0;
  }

  /* To avoid problems with incomplete data bases from crashed
     calculations, delete the data base if we can ... but DON'T
     do this until we have checked that it is not already open */

  if ((exists && (flags & O_TRUNC)) ||
      (strcmp(mode, "unknown") == 0 && file_size(filename) <= 0))
    (void) unlink(filename);

  /* Open the physical data base */

#if USE_DB  
  openinfo.bsize = 1024*4;
  openinfo.cachesize = 128*1024;
  openinfo.ffactor = 16;
  openinfo.hash = 0;
  openinfo.lorder = 0;
  openinfo.nelem = 1024;

  if (!(rtdb[new].db = 
	dbopen(filename, flags, 0660, DB_HASH, &openinfo))) {
    (void) fprintf(stderr, "rtdb_seq_open: db failed to open file %s\n",
		   filename);
    return 0;
  }
#endif

#if USE_HDBM
#ifdef DB_INMEM
  if (!hdbm_open(filename, 0, &rtdb[new].db)) {
#else
  if (!hdbm_open(filename, 1, &rtdb[new].db)) {
#endif
    (void) fprintf(stderr, "rtdb_seq_open: hdbm failed to open file %s\n",
		   filename);
    return 0;
  }
#endif

  /* Shove info into the RTDB structure and have at it */

  rtdb[new].active = 1;
  rtdb[new].filename = strdup(filename);
  rtdb[new].scratch = (strcmp(mode, "scratch") == 0);

  *handle = new;

  return 1;
}

static int check_handle(const int handle)
{
  if (handle < 0 || handle >= MAX_RTDB) {
    (void) fprintf(stderr, "check_handle: handle (%d) is out of range\n",
		   handle);
    return 0;
  }
  
  if (!rtdb[handle].active) {
    (void) fprintf(stderr, "check_handle: handle (%d) is inactive\n", handle);
    return 0;
  }

  return 1;
}

int rtdb_seq_first(const int handle, const int namelen, char *name)
/*
  Return the name of the first (user inserted) entry in the data base.
  The order is effectively random.

  handle  = handle to RTDB
  namelen = size of user provided buffer name
  name    = name of entry is returned in this buffer
*/
{
#if USE_DB  
  DB *db;
  int flags = R_FIRST;
#endif
#if USE_HDBM
  hdbm db;
#endif
  int status;
  DBT key;

#ifdef USE_DB
  DBT value;
#endif

  if (!check_handle(handle)) {
    (void) fprintf(stderr, "rdtb_first: handle (%d) is invalid\n", handle);
    return 0;
  }

  db = rtdb[handle].db;

#ifdef USE_DB
  while ((status = db->seq(db, &key, &value, flags)) == 0) {
    flags = R_NEXT;
    if (strncmp(info_header, DATA(key), strlen(info_header)) != 0)
      break;
  }
  status = !status;
#endif
#ifdef USE_HDBM
  for (status=hdbm_first_key(db, &key);
       status;
       status=hdbm_next_key(db, &key)) {
    if (strncmp(info_header, DATA(key), strlen(info_header)) != 0)
      break;
    datum_free(key);
  }
#endif

  if (status) {
    if (namelen >= SIZE(key)) {
      strncpy(name, (char *) DATA(key), namelen);
      name[namelen-1] = 0;
    }
    else {
      (void) fprintf(stderr, 
		     "rtdb_seq_first: name is too small, need=%u got=%d\n",
		     (int) SIZE(key), (int) namelen);
      status = 0;
      }
#ifdef USE_HDBM
    datum_free(key);
#endif
  }
  return status;
}

int rtdb_seq_delete(const int handle, const char *name)
/*
  Delete the entry from the database.
  Return
        1 if key was present and successfully deleted

	0 if key was not present, or if an error occured

  handle  = handle to RTDB
  name    = name of entry to delete
*/
{
#ifdef USE_DB
  DB *db = rtdb[handle].db;
#endif
#ifdef USE_HDBM
  hdbm db = rtdb[handle].db;
#endif
  DBT key, info_key;
  int status;
  char info_buf[256];

#ifdef USE_DB  
  int flags = 0;
#endif

  if (!check_handle(handle)) {
    (void) fprintf(stderr, "rdtb_delete: handle (%d) is invalid\n", handle);
    return 0;
  }

  if (strlen(name)+sizeof(info_header)+1 > sizeof(info_buf)) {
    (void) fprintf(stderr,
		   "rtdb_delete: info entry buffer needs to be %d\n",
		   (int) (strlen(name)+sizeof(info_header)+1));
    return 0;
  }
  
  (void) sprintf(info_buf, "%s%s", info_header, name);

  WRAP(info_buf, strlen(info_buf)+1, &info_key);
  WRAP(name, strlen(name)+1, &key);

#ifdef USE_DB  
  status = db->del(db, &key, flags);
  if (status == 0)
    status = db->del(db, &info_key, flags);

#endif
#ifdef USE_HDBM
  status = !hdbm_delete(db, key);
  if (status == 0)
    status = !hdbm_delete(db, info_key);
#endif

  if (status == 0)
    return 1;
  else if (status == 1)
    return 0;
  else {
    (void) fprintf(stderr, "rtdb_seq_delete: error deleting %s\n", name);
    return 0;
  }
}
    
int rtdb_seq_next(const int handle, const int namelen, char *name)
/*
  Return the name of the next (user inserted) entry in the data base.
  The order is effectively random.

  handle  = handle to RTDB
  namelen = size of user provided buffer name
  name    = name of entry is returned in this buffer
*/
{
#ifdef USE_DB
  DB *db = rtdb[handle].db;
  int flags = R_NEXT;
#endif
#ifdef USE_HDBM
  hdbm db = rtdb[handle].db;
#endif
  int status;
  DBT key;

#ifdef USE_DB
  DBT value;
#endif

  if (!check_handle(handle)) {
    (void) fprintf(stderr, "rdtb_close: handle (%d) is invalid\n", handle);
    return 0;
  }

#ifdef USE_DB
  while ((status = db->seq(db, &key, &value, flags)) == 0) { 
    flags = R_NEXT;		/* } For emacs */
    if (strncmp(info_header, DATA(key), strlen(info_header)) != 0)
      break;
  }
  status = !status;
#endif
#ifdef USE_HDBM
  while ((status = hdbm_next_key(db, &key))) {
    if (strncmp(info_header, DATA(key), strlen(info_header)) != 0)
      break;
    datum_free(key);
  }
#endif

  if (status) {
    if (namelen >= SIZE(key)) {
      strncpy(name, (char *) DATA(key), namelen);
    }
    else {
      (void) fprintf(stderr, 
		     "rtdb_seq_next: name too small, need=%u, got=%d\n",
		     (int) SIZE(key), (int) namelen);
      status = 0;
    }
#ifdef USE_HDBM
    datum_free(key);
#endif
  }

  return status;
}


#include "ga.h"
int rtdb_seq_close(const int handle, const char *mode)
/*
  Close the data base

  handle  = handle to RTDB
  mode    = 'keep'    Preserve the data base file to enable restart
            'delete'  Delete the data base file freeing all resources

  mode is overridden by opening the data base with mode='scratch'
  in which instance it is always deleted upon closing
*/
{
#ifdef USE_DB
  DB *db = rtdb[handle].db;
#endif
#ifdef USE_HDBM
  hdbm db = rtdb[handle].db;
#endif

  if (!check_handle(handle)) {
    (void) fprintf(stderr, "rdtb_close: handle (%d) is invalid\n", handle);
    return 0;
  }

  /* Mark as inactive even if we trap any other errors */

  rtdb[handle].active = 0;

#ifdef USE_DB
  /* Flush changes to the DB if required */

  if (!(strcmp(mode,"delete") == 0 || rtdb[handle].scratch)) {

    if (db->sync(db, (u_int) 0)) {
      (void) fprintf(stderr, "rdtb_close: db sync(%s) returned error\n", 
		     rtdb[handle].filename);
      return 0;
    }
  }  

  /* Actually close the data base */
  
  if (db->close(db)) {
    (void) fprintf(stderr, "rdtb_close: db close(%s) returned error\n", 
		   rtdb[handle].filename);
    return 0;
  }
#endif
#ifdef USE_HDBM
  if (!hdbm_close(db)) {
    (void) fprintf(stderr, "rdtb_close: hdbm close(%s) returned error\n", 
		   rtdb[handle].filename);
    return 0;
  }
#endif

  /* Delete it if required */
    
  if (strcmp(mode,"delete") == 0 || rtdb[handle].scratch) {
    if (unlink(rtdb[handle].filename) != 0) {
      (void) fprintf(stderr, "rdtb_close: error in deleting %s\n",
		     rtdb[handle].filename);
      return 0;
    }
  }
#ifdef USE_HDBM
  else {			/* Compress out dead space */
    if (!hdbm_file_compress(rtdb[handle].filename)) {
      (void) fprintf(stderr, "rdtb_close: hdbm compress(%s) returned error\n", 
		     rtdb[handle].filename);
      return 0;
    }
  }
#endif

  return 1;
}

static void get_time(char buf[26])
{
  time_t t = time((time_t *) 0);
  char *tmp = ctime(&t);

  (void) memcpy(buf, tmp, 26);
}

int rtdb_seq_print(const int handle, const int print_values)
/*
  Print the contents of the data base to stdout

  handle  = handle to RTDB
  print_values = boolean flag ... if true values as well as
                 keys are printed out.
*/
{
  char name[128];
  int status;
  int len;

  if (!check_handle(handle)) {
    (void) fprintf(stderr, "rtdb_seq_print: handle (%d) is invalid\n", handle);
    return 0;
  }

  printf("\n Contents of RTDB %s\n ", rtdb[handle].filename);
  len = strlen(rtdb[handle].filename) + 17;
  while(len--)
    printf("-");
  printf("\n\n");
  printf(" Entry                                   Type[nelem]  "
	 "         Date\n");
  printf(" ---------------------------  ----------------------  "
	 "------------------------\n");

  status = rtdb_seq_first(handle, sizeof(name), name);

  while (status) {
    char date[26];
    int nelem;
    int ma_type;

    if (!rtdb_seq_get_info(handle, name, &ma_type, &nelem, date)) {
      (void) fprintf(stderr, "rtdb_seq_print: get_info failed on %s\n", name);
    }
    else {			/* Fortran is better at this ! */
      printf(" %s", name);
      len = 29 - strlen(name);
      while (len-- > 0)
	printf(" ");
      printf("%15s[%d]", ma_typename(ma_type), nelem);
      if (nelem < 10)
	printf(" ");
      if (nelem < 100)
	printf(" ");
      if (nelem < 1000)
	printf(" ");
      if (nelem < 10000)
	printf(" ");
      if (nelem < 100000)
	printf(" ");
      printf(" %s\n", date);

      if (print_values) {
	int ma_handle;

	if (!rtdb_seq_ma_get(handle, name, &ma_type, &nelem, &ma_handle)) {
	  (void) fprintf(stderr, "rtdb_seq_print: ma_get failed on %s\n", name);
	}
	else {
	  void *data;
	  if (MA_get_pointer(ma_handle, &data))
	    ma_print(stdout, ma_type, nelem, data);
	  else
	    (void) fprintf(stderr, "rtdb_seq_print: MA_get_pt failed, handle=%d\n",
			   handle);
	  (void) MA_free_heap(ma_handle);
	}
      }
    }
    status = rtdb_seq_next(handle, sizeof(name), name);
  }
  printf("\n");

  return 1;
}

static int rtdb_seq_put_info(const int handle,
			 const char *name, const int ma_type, const int nelem)
/*
  Insert info about an entry into the data base replacing previous entry

  handle   = handle to RTDB
  name     = entry name (null terminated character string)
  ma_type  = MA type of the entry
  nelem    = no. of elements of the given type
*/  
{
  struct info_struct info;
  char info_buf[256];
#ifdef USE_DB
  DB *db = rtdb[handle].db;
#endif
#ifdef USE_HDBM
  hdbm db = rtdb[handle].db;
#endif
  DBT key, value;
  int status;

  if (!check_handle(handle)) {
    (void) fprintf(stderr, "rtdb_seq_put_info: handle (%d) is invalid\n", handle);
    return 0;
  }

  (void) memset((void *) &info, (char) 0, 
		(size_t) (sizeof(struct info_struct)));

  info.ma_type = ma_type;
  info.nelem   = nelem;
  get_time(info.date);

  if (strlen(name)+sizeof(info_header)+1 > sizeof(info_buf)) {
    (void) fprintf(stderr,
		   "rtdb_seq_put_info: info entry buffer needs to be %d\n",
		   (int) (strlen(name)+sizeof(info_header)+1));
    return 0;
  }
  
  (void) sprintf(info_buf, "%s%s", info_header, name);

  WRAP(info_buf, strlen(info_buf)+1, &key);
  WRAP(&info, sizeof(info), &value);
  
#ifdef USE_DB
  status = !(db->put(db, &key, &value, (u_int) 0));
#endif
#ifdef USE_HDBM
  status = hdbm_replace(db, key, value);
#endif

  if (!status) {
    (void) fprintf(stderr, 
		   "rtdb_seq_put_info: put failed for \"%s\" in %s\n",
		   name, rtdb[handle].filename);
    return 0;
  }

  return 1;
}

int rtdb_seq_get_info(const int handle,
		  const char *name, int *ma_type, int *nelem, char date[26])
/*
  Get info about an entry from the data base

  handle   = handle to RTDB
  name     = entry name (null terminated character string)
  ma_type  = returns MA type of the entry
  nelem    = returns no. of elements of the given type
  date     = returns date of insertion (null terminated character string)
*/  
{
  struct info_struct info;
  char info_buf[256];
#ifdef USE_DB
  DB *db = rtdb[handle].db;
#endif
#ifdef USE_HDBM
  hdbm db = rtdb[handle].db;
#endif
  DBT key, value;
  int status;

  if (!check_handle(handle)) {
    (void) fprintf(stderr, "rtdb_seq_get_info: handle (%d) is invalid\n", 
		   handle);
    return 0;
  }

  (void) memset((void *) &info, (char) 0, 
		(size_t) (sizeof(struct info_struct)));


  if (strlen(name)+sizeof(info_header)+1 > sizeof(info_buf)) {
    (void) fprintf(stderr,
		   "rtdb_seq_get_info: info entry buffer needs to be %d\n",
		   (int) (strlen(name)+sizeof(info_header)+1));
    return 0;
  }
  
  (void) sprintf(info_buf, "%s%s", info_header, name);

  WRAP(info_buf, strlen(info_buf)+1, &key);

#ifdef USE_DB
  status = (db->get(db, &key, &value, (u_int) 0));
  if (status == -1) {
    (void) fprintf(stderr, 
		   "rtdb_seq_get_info: db get failed for \"%s\" in %s\n",
		   name, rtdb[handle].filename);
  }
  status = !status;
#endif
#ifdef USE_HDBM
  status = hdbm_read(db, key, &value);
#endif

  if (!status) {
      /* Entry not found ... quietly return error so that failed
	 probes are not excessively verbose */
    return 0;
  }

  if (SIZE(value) != sizeof(info)) {
    (void) fprintf(stderr, 
		   "rtdb_seq_get_info: size mismatch : info=%d, db=%d\n",
		   (int) sizeof(info), (int) SIZE(value));
    return 0;
  }

  (void) memcpy(&info, DATA(value), SIZE(value));
#ifdef USE_HDBM
  datum_free(value);
#endif

  *ma_type = info.ma_type;
  *nelem   = info.nelem;
  if (info.date[24] == '\n') info.date[24] = ' ';
  (void) memcpy(date, info.date, 26);

  return 1;
}

int rtdb_seq_put(const int handle, const char *name, const int ma_type,
	     const int nelem, const void *array)
/*
  Insert an entry into the data base replacing previous entry

  handle   = handle to RTDB
  name     = entry name (null terminated character string)
  ma_type  = MA type of the entry
  nelem    = no. of elements of the given type
  array    = data to be inserted
*/  
{
#ifdef USE_DB
  DB *db = rtdb[handle].db;
#endif
#ifdef USE_HDBM
  hdbm db = rtdb[handle].db;
#endif
  DBT key, value;

  int status;

  if (!check_handle(handle)) {
    (void) fprintf(stderr, "rtdb_seq_put: handle (%d) is invalid\n", handle);
    return 0;
  }

  /* Enter the data into the data base */

  WRAP(name, strlen(name)+1,&key);
  WRAP(array, MA_sizeof(ma_type, nelem, MT_CHAR),&value);

#ifdef USE_DB
  status = !db->put(db, &key, &value, (u_int) 0);
#endif
#ifdef USE_HDBM
  status = hdbm_replace(db, key, value);
#endif

  if (!status) {
    (void) fprintf(stderr, "rtdb_seq_put: put failed for \"%s\" in %s\n",
		   name, rtdb[handle].filename);
    return 0;
  }

  /* Enter the info into the data base as "!rtdb!<name>" */

  if (!rtdb_seq_put_info(handle, name, ma_type, nelem))
    return 0;

  /* This operations flushes ALL writes immediately to disk
     ... it may be slow, in which case comment the next statement
     out and be sure to explicitly close the databse */

#ifdef USE_DB
  status = (db->sync(db, (u_int) 0) == 0);
#endif
#ifdef USE_HDBM
#ifndef DB_INMEM
  status = hdbm_file_flush(db);  
#endif
#endif

  if (!status) {
    (void) fprintf(stderr, "rtdb_seq_put: flush failed for \"%s\" in %s\n",
		   name, rtdb[handle].filename);
    return 0;
  }

  return 1;
}

int rtdb_seq_get(const int handle, const char *name, const int ma_type,
	     const int nelem, void *array)
/*
  Get an entry from the data base

  handle   = handle to RTDB
  name     = entry name (null terminated character string)
  ma_type  = MA type of the entry which must match entry type
  nelem    = size of array in units of ma_type
  array    = user provided buffer that returns data
*/
{
#ifdef USE_DB
  DB *db = rtdb[handle].db;
#endif
#ifdef USE_HDBM
  hdbm db = rtdb[handle].db;
#endif

  DBT key, value;
  int status, ma_type_from_db, nelem_from_db;
  char date[26];
  int nelem_actual;

  if (!check_handle(handle)) {
    (void) fprintf(stderr, "rtdb_seq_get: handle (%d) is invalid\n", handle);
    return 0;
  }

  /* Retrieve the info from the data base to check types etc. */

  if (!rtdb_seq_get_info(handle, name, &ma_type_from_db, &nelem_from_db, date)) {

    /* In production will want to have this handled quietly  ... be
       verbose for debug only 

    (void) fprintf(stderr, "rtdb_seq_get: get info failed for \"%s\" in %s\n",
		   name, rtdb[handle].filename); */
    return 0;
  }

  if (ma_type_from_db != ma_type) {
    (void) fprintf(stderr, 
		   "rtdb_seq_get: type mismatch \"%s\" in %s: arg=%d, db=%d\n",
		   name, rtdb[handle].filename, ma_type, ma_type_from_db);
    return 0;
  }

  /* Get the data from the data base */

  WRAP(name, strlen(name)+1, &key);
#ifdef USE_DB
  /* Return can be 1 (not present) or -1 (failed).  Both are an
     error here as the get_info above succeeded */
  status = (db->get(db, &key, &value, (u_int) 0) == 0);
#endif
#ifdef USE_HDBM
  status = hdbm_read(db, key, &value);
#endif

  if (!status) {
    (void) fprintf(stderr, "rtdb_seq_get: get failed for \"%s\" in %s\n",
		   name, rtdb[handle].filename);
    return 0;
  }

  /* Now check that user's buffer is big enough */

  nelem_actual = MA_sizeof(MT_CHAR, SIZE(value), ma_type);

  if (nelem_actual != nelem_from_db) {
    (void) fprintf(stderr, 
		   "rtdb_seq_get: size error \"%s\" in %s: info=%d, db=%d\n",
		   name, rtdb[handle].filename, nelem_from_db,
		   nelem_actual);
    return 0;
  }

  if (nelem_actual > nelem) {
    (void) fprintf(stderr, "rtdb_seq_get: \"%s\" is %d long, buffer is %d\n",
		   name, nelem_actual, nelem);
    return 0;
  }

  (void) memcpy(array, DATA(value), SIZE(value));
#ifdef USE_HDBM
  datum_free(value);
#endif

  return 1;
}

int rtdb_seq_ma_get(const int handle, const char *name, int *ma_type,
		int *nelem, int *ma_handle)
/*
  Get an entry from the data base returning an MA handle

  handle   = handle to RTDB
  name     = entry name (null terminated character string)
  ma_type  = returns MA type of the entry
  nelem    = returns no. of elements of type ma_type in data
  ma_handle= returns MA handle to data
*/
{
  char date[26];
  void *ma_data;
  Integer ma_handle_buf;
  Integer ma_type_buf;
  Integer nelem_buf;

  if (!check_handle(handle)) {
    (void) fprintf(stderr, "rtdb_seq_ma_get: handle (%d) invalid\n", handle);
    return 0;
  }

  /* Retrieve the type info from the data base */

  if (!rtdb_seq_get_info(handle, name, ma_type, nelem, date)) {
    /*
    (void) fprintf(stderr, "rtdb_seq_ma_get: get info failed \"%s\" in %s\n",
		   name, rtdb[handle].filename);
                   */
    return 0;
  }

  /* Allocate MA pointers */

  ma_type_buf = (Integer) *ma_type;
  nelem_buf   = (Integer) *nelem;

  if (!MA_allocate_heap(ma_type_buf, nelem_buf, name, &ma_handle_buf)) {
    (void) fprintf(stderr, "rtdb_seq_ma_get: MA_allocate_heap : nelem=%d\n",
		   *nelem);
    return 0;
  }
  *ma_handle = ma_handle_buf;
  
  if (!MA_get_pointer(ma_handle_buf, &ma_data)) {
    (void) fprintf(stderr, "rtdb_seq_ma_get: MA_get_pointer failed\n");
    return 0;
  }

  if (!rtdb_seq_get(handle, name, *ma_type, *nelem, ma_data)) {
    (void) fprintf(stderr, "rtdb_seq_ma_get: rtdb_seq_get failed %s\n", name);
    (void) MA_free_heap(ma_handle_buf);
    return 0;
  }

  return 1;
}

int rtdb_seq_copy(const int handle, const char *suffix)
/*
  Copy the data base

  handle  = handle to RTDB
  suffix    = new file will be called rtdb_name.suffix
*/
{

#ifdef USE_HDBM
    if (!hdbm_file_copy(rtdb[handle].filename, suffix)) {
	(void) fprintf(stderr,
		       "rtdb_seq_copy: copy from %s failed\n", suffix);
      return 0;
    }
#else
fixme
#endif
  return 1;
}
int rtdb_seq_getfname(const int handle,
		  char* fname)
/*
  Get rtdb file name

  handle   = handle to RTDB
  fname    = returns rtdb file name (null terminated character string)
*/  
{

  if (!check_handle(handle)) {
    (void) fprintf(stderr, "rtdb_seq_get_info: handle (%d) is invalid\n", 
		   handle);
    return 0;
  }

  (void) memcpy(fname, rtdb[handle].filename, strlen(rtdb[handle].filename)+1);
/*  (void) memcpy(fname, rtdb[handle].filename, 36); */

  return 1;
}

