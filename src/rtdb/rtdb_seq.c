#include <stdlib.h>
#include <sys/types.h>
#include <stdio.h>
#include <db.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#ifndef IPSC
#include <sys/time.h>
#endif
/*#include <limits.h>*/
#include <sys/stat.h>
#include <fcntl.h>
#include "rtdb.h"
#include "macdecls.h"
#include "misc.h"

#include "sndrcv.h"

extern char *strdup(const char *);
extern void *malloc(size_t);
extern void free(void *);

#define MAX_RTDB 5

#define FORTRAN_TRUE  1L
#define FORTRAN_FALSE 0L

static struct {			/* Keep track of active RTDBs */
  int active;
  char *filename;
  int scratch;
  DB *db;
} rtdb[MAX_RTDB];

struct info_struct{		/* Matching info entry for each datum */
  int ma_type;
  int nelem;
  char date[26];
};

static char info_header[] = "!rtdb!"; /* Prefix for info entries */

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
  case MT_BASE + 0:	/* char */
  case MT_BASE + 9:	/* Fortran byte */

    (void) fprintf(file, "%.*s\n", nelem, (char *) p);
    break;

  case MT_BASE + 1:	/* int */
    for (nprint=i=0; i<nelem; i++) {
      nprint += fprintf(file, "%d ", ((int *) p)[i]);
      if (nprint >= 72) {
	(void) fprintf(file, "\n");
	nprint = 0;
      }
    }
    if (nprint > 0) (void) fprintf(file, "\n");
    break;

  case MT_BASE + 10:	/* Fortran integer ... not equivalent on KSR */
    for (nprint=i=0; i<nelem; i++) {
      nprint += fprintf(file, "%d ", ((Integer *) p)[i]);
      if (nprint >= 72) {
	(void) fprintf(file, "\n");
	nprint = 0;
      }
    }
    if (nprint > 0) (void) fprintf(file, "\n");
    break;


  case MT_BASE + 2:	/* long int */

    for (nprint=i=0; i<nelem; i++) {
      nprint += fprintf(file, "%ld ", ((long *) p)[i]);
      if (nprint >= 72) {
	(void) fprintf(file, "\n");
	nprint = 0;
      }
    }
    if (nprint > 0) (void) fprintf(file, "\n");
    break;

  case MT_BASE + 3:	/* float */
  case MT_BASE + 12:	/* Fortran real */

    for (nprint=i=0; i<nelem; i++) {
      nprint += fprintf(file, "%g ", ((float *) p)[i]);
      if (nprint >= 72) {
	(void) fprintf(file, "\n");
	nprint = 0;
      }
    }
    if (nprint > 0) (void) fprintf(file, "\n");
    break;

  case MT_BASE + 4:	/* double */
  case MT_BASE + 13:	/* Fortran double precision */

    for (nprint=i=0; i<nelem; i++) {
      nprint += fprintf(file, "%g ", ((double *) p)[i]);
      if (nprint >= 72) {
	(void) fprintf(file, "\n");
	nprint = 0;
      }
    }
    if (nprint > 0) (void) fprintf(file, "\n");
    break;

  case MT_BASE + 6:	/* single precision complex */
  case MT_BASE + 14:	/* Fortran single precision complex */

    for (nprint=i=0; i<nelem; i++) {
      nprint += fprintf(file, "(%g,%g) ", ((float *) p)[2*i], 
			((float *) p)[2*i+1]);
      if (nprint >= 72) {
	(void) fprintf(file, "\n");
	nprint = 0;
      }
    }
    if (nprint > 0) (void) fprintf(file, "\n");
    break;

  case MT_BASE + 7:	/* double precision complex */
  case MT_BASE + 15:	/* Fortran double precision complex */

    for (nprint=i=0; i<nelem; i++) {
      nprint += fprintf(file, "(%g,%g) ", ((double *) p)[2*i], 
			((double *) p)[2*i+1]);
      if (nprint >= 72) {
	(void) fprintf(file, "\n");
	nprint = 0;
      }
    }
    if (nprint > 0) (void) fprintf(file, "\n");
    break;

  case MT_BASE + 11:	/* Fortran logical */
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


  case MT_BASE + 8:	/* long double precision complex */
  case MT_BASE + 5:	/* long double */

    (void) fprintf(file, " !! printing long double not supported !!\n");
    break;

  default:

    (void) fprintf(file, " !! %d = unknown data type\n", ma_type);
  }
}

const char *ma_typename(const int ma_type)
{
  switch (ma_type) {
  case MT_BASE + 0:	/* char */
  case MT_BASE + 9:	/* Fortran byte */
    return "char"; break;

  case MT_BASE + 1:	/* int */
  case MT_BASE + 10:	/* Fortran integer ... not equivalent on KSR */

    return "int"; break;

  case MT_BASE + 2:	/* long int */

    return "long"; break;

  case MT_BASE + 3:	/* float */
  case MT_BASE + 12:	/* Fortran real */

    return "float"; break;

  case MT_BASE + 4:	/* double */
  case MT_BASE + 13:	/* Fortran double precision */

    return "double"; break;

  case MT_BASE + 6:	/* single precision complex */
  case MT_BASE + 14:	/* Fortran single precision complex */

    return "complex float";

  case MT_BASE + 7:	/* double precision complex */
  case MT_BASE + 15:	/* Fortran double precision complex */

    return "complex double"; break;

  case MT_BASE + 11:	/* Fortran logical */

    return "logical"; break;

  case MT_BASE + 8:	/* long double precision complex */
  case MT_BASE + 5:	/* long double */
  default:

    return "invalid"; break;
  }
}


static DBT make_DBT(const void *data, size_t size) 
/*
  Make a DBT by duplicating the data described
  by the arguments.
*/
{
  DBT tmp;

  tmp.size = 0; tmp.data = 0;
  
  if (size && data) {
    ALLOCATE(tmp.data, size, char, "make_DBT: tmp.data (%d)");
    memcpy(tmp.data, (const void *) data, (int) size);
    tmp.size = size;
  }

  return tmp;
}

static void print_hex_DBT(const DBT *d)
/*
  Print summary of DBT to stdout in hex
*/
{
  long len = MIN(8, d->size);
  long i;
  unsigned char *ptr = (unsigned char *) d->data;

  printf("(%d, 0x%p, 0x",(int) d->size, d->data);
  if (d->data) {
    for (i=0; i<len; i++) 
      printf("%2.2x", (unsigned) ptr[i]);
    if (len < d->size)
      printf(" ...");
  }
  printf(")");
  fflush(stdout);
}

static void print_string_DBT(const DBT *d)
/*
  Prints DBT to stdout assuming a null terminated string.
*/
{
  long len = MIN(8, d->size-1);

  (void) printf("(%ld, 0x%lx, \"",d->size, (unsigned long) d->data);
  if (d->data) {
    (void) printf("%.*s", (int) len, (char *) d->data);
    if (len < (d->size-1))
      (void) printf("...");
  }
  (void) printf("\")");
  fflush(stdout);
}

static DBT duplicate_DBT(const DBT *d)
/*
  Duplicate the DBT d, ALLOCATing space as required.
*/
{
  void *data = d->data;
  size_t size = d->size;
  DBT tmp;

  tmp.data = 0; tmp.size = 0;
  
  if (size && data) {
    ALLOCATE(tmp.data, size, char, "duplicate_DBT: data (%d)");
    memcpy(tmp.data, data, size);
    tmp.size = size;
  }

  return tmp;
}

static void free_DBT(DBT *d)
/*
  Free the memory pointed to by d.data, assumed to have
  got from ALLOCATE.
*/
{
  if (d->data)
    FREE(d->data);
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
  int flags = O_RDWR | O_CREAT;
  int exists = access(filename, R_OK | W_OK) == 0;
  int new;
  HASHINFO openinfo;

  /* See if the data base is already open ... if so return and complain */

  for (new=0; new<MAX_RTDB; new++)
    if (rtdb[new].active)
      if (strcmp(filename, rtdb[new].filename) == 0) {
	(void) fprintf(stderr, "rtdb_seq_open: %s is already open\n", filename);
	return 0;
      }

  /* Figure out file access modes */

  if (strcmp(mode, "new") == 0) {
    if (exists) {
      (void) fprintf(stderr, "rtdb_seq_open: %s exists but is being opened new\n",
		     filename);
      return 0;
    }
  }
  else if (strcmp(mode, "old") == 0) {
    if (!exists) {
      (void) fprintf(stderr, 
		     "rtdb_seq_open: %s does not exist but is being opened old\n",
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

  openinfo.bsize = 1024*4;
#if defined(CRAY)
  openinfo.cachesize = 4*128*1024;
#else
  openinfo.cachesize = 128*1024;
#endif
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
  DB *db;
  DBT key, value;
  int flags = R_FIRST, status;

  if (!check_handle(handle)) {
    (void) fprintf(stderr, "rdtb_first: handle (%d) is invalid\n", handle);
    return 0;
  }

  db = rtdb[handle].db;

  while ((status = db->seq(db, &key, &value, flags)) == 0) {

    flags = R_NEXT;

    if (strncmp(info_header, key.data, strlen(info_header)) != 0)
      break;
  }

  if (status == 0) {
    if (namelen >= key.size) {
      strncpy(name, (char *) key.data, namelen);
      name[namelen-1] = 0;
      return 1;
    }
    else {
      (void) fprintf(stderr, "rtdb_seq_first: name is too small, need=%u got=%d\n",
		     (int) key.size, (int) namelen);
      return 0;
    }
  }
  else
    return 0;
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
  DB *db;
  DBT key;
  int flags = 0, status;

  if (!check_handle(handle)) {
    (void) fprintf(stderr, "rdtb_delete: handle (%d) is invalid\n", handle);
    return 0;
  }

  key = wrap_DBT(name, strlen(name)+1);

  db = rtdb[handle].db;

  status = db->del(db, &key, flags);

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
  DB *db;
  DBT key, value;
  int flags = R_NEXT, status;

  if (!check_handle(handle)) {
    (void) fprintf(stderr, "rdtb_close: handle (%d) is invalid\n", handle);
    return 0;
  }

  db = rtdb[handle].db;

  while ((status = db->seq(db, &key, &value, flags)) == 0) {
    flags = R_NEXT;

    if (strncmp(info_header, key.data, strlen(info_header)) != 0)
      break;
  }

  if (status == 0) {
    if (namelen >= key.size) {
      strncpy(name, (char *) key.data, namelen);
      return 1;
    }
    else {
      (void) fprintf(stderr, "rtdb_seq_next: name too small, need=%u, got=%d\n",
		     (int) key.size, (int) namelen);
      return 0;
    }
  }
  else
    return 0;
}


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
  DB *db;

  if (!check_handle(handle)) {
    (void) fprintf(stderr, "rdtb_close: handle (%d) is invalid\n", handle);
    return 0;
  }

  db = rtdb[handle].db;
  
  /* Mark as inactive even if we trap any other errors */

  rtdb[handle].active = 0;

  /* Flush changes to the DB if required */

  if (!(strcmp(mode,"delete") == 0 || rtdb[handle].scratch))
    if (NODEID_() == 0)
      if (db->sync(db, (u_int) 0)) {
	(void) fprintf(stderr, "rdtb_close: db sync(%s) returned error\n", 
		       rtdb[handle].filename);
	return 0;
      }
  
  /* Actually close the data base */
  
  if (db->close(db)) {
    (void) fprintf(stderr, "rdtb_close: db close(%s) returned error\n", 
		   rtdb[handle].filename);
    return 0;
  }

  /* Delete it if required */
    
  if (strcmp(mode,"delete") == 0 || rtdb[handle].scratch) {
    if (NODEID_() == 0) {
      if (unlink(rtdb[handle].filename) != 0) {
	(void) fprintf(stderr, "rdtb_close: error in deleting %s\n",
		       rtdb[handle].filename);
	return 0;
      }
    }
  }

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
      printf(" %s", date);

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
  DB *db = rtdb[handle].db;
  DBT key, value;

  if (!check_handle(handle)) {
    (void) fprintf(stderr, "rtdb_seq_put_info: handle (%d) is invalid\n", handle);
    return 0;
  }

  info.ma_type = ma_type;
  info.nelem   = nelem;
  get_time(info.date);

  if (strlen(name)+sizeof(info_header)+1 > sizeof(info_buf)) {
    (void) fprintf(stderr,"rtdb_seq_put_info: info entry buffer needs to be %d\n",
		   (int) (strlen(name)+sizeof(info_header)+1));
    return 0;
  }
  
  (void) sprintf(info_buf, "%s%s", info_header, name);

  key = wrap_DBT(info_buf, strlen(info_buf)+1);
  value = wrap_DBT(&info, sizeof(info));

  if (db->put(db, &key, &value, (u_int) 0)) {
    (void) fprintf(stderr, "rtdb_seq_put_info: db put failed for \"%s\" in %s\n",
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
  DB *db = rtdb[handle].db;
  DBT key, value;
  int status;

  if (!check_handle(handle)) {
    (void) fprintf(stderr, "rtdb_seq_get_info: handle (%d) is invalid\n", handle);
    return 0;
  }

  if (strlen(name)+sizeof(info_header)+1 > sizeof(info_buf)) {
    (void) fprintf(stderr,"rtdb_seq_get_info: info entry buffer needs to be %d\n",
		   (int) (strlen(name)+sizeof(info_header)+1));
    return 0;
  }
  
  (void) sprintf(info_buf, "%s%s", info_header, name);

  key = wrap_DBT(info_buf, strlen(info_buf)+1);

  if (status = db->get(db, &key, &value, (u_int) 0)) {

    if (status == -1) {
      (void) fprintf(stderr, "rtdb_seq_get_info: db get failed for \"%s\" in %s\n",
		     name, rtdb[handle].filename);
    }
    else {
      /* Entry not found ... quietly return error so that failed
	 probes are not excessively verbose */
    }
    return 0;
  }

  if (value.size != sizeof(info)) {
    (void) fprintf(stderr, "rtdb_seq_get_info: size mismatch : info=%d, db=%d\n",
		   (int) sizeof(info), (int) value.size);
    return 0;
  }

  (void) memcpy(&info, value.data, value.size);

  *ma_type = info.ma_type;
  *nelem   = info.nelem;
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
  DBT key, value;
  DB *db;
  int status;

  if (!check_handle(handle)) {
    (void) fprintf(stderr, "rtdb_seq_put: handle (%d) is invalid\n", handle);
    return 0;
  }

  /* Enter the data into the data base */

  key = wrap_DBT(name, strlen(name)+1);
  value = wrap_DBT(array, MA_sizeof(ma_type, nelem, MT_CHAR));

  db = rtdb[handle].db;

  if (db->put(db, &key, &value, (u_int) 0)) {
    (void) fprintf(stderr, "rtdb_seq_put: db put failed for \"%s\" in %s\n",
		   name, rtdb[handle].filename);
    return 0;
  }

  /* Enter the info into the data base as "!rtdb!<name>" */

  status = rtdb_seq_put_info(handle, name, ma_type, nelem);

  /* This operations flushes ALL writes immediately to disk
     ... it may be slow, in which case comment the next statement
     out and be sure to explicitly close the databse */

  status = status && (db->sync(db, (u_int) 0) == 0);

  return status;
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
  DBT key, value;
  int status, ma_type_from_db, nelem_from_db;
  DB *db;
  char date[26];
  int nelem_actual;

  if (!check_handle(handle)) {
    (void) fprintf(stderr, "rtdb_seq_get: handle (%d) is invalid\n", handle);
    return 0;
  }

  db = rtdb[handle].db;

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
		   "rtdb_seq_get: type mismatch for \"%s\" in %s: arg=%d, db=%d\n",
		   name, rtdb[handle].filename, ma_type, ma_type_from_db);
    return 0;
  }

  /* Get the data from the data base */

  key = wrap_DBT(name, strlen(name)+1);

  if ((status = db->get(db, &key, &value, (u_int) 0))) {
    /* Status can be 1 (not present) or -1 (failed).  Both are an
       error here as the get_info above succeeded */
    
    (void) fprintf(stderr, "rtdb_seq_get: db get failed for \"%s\" in %s\n",
		   name, rtdb[handle].filename);
    return 0;
  }

  /* Now check that users buffer is big enough */

  nelem_actual = MA_sizeof(MT_CHAR, value.size, ma_type);

  if (nelem_actual != nelem_from_db) {
    (void) fprintf(stderr, 
		   "rtdb_seq_get: size error for \"%s\" in %s: info=%d, db=%d\n",
		   name, rtdb[handle].filename, nelem_from_db,
		   nelem_actual);
    return 0;
  }

  if (nelem_actual > nelem) {
    (void) fprintf(stderr, "rtdb_seq_get: \"%s\" is %d long, buffer is only %d\n",
		   name, nelem_actual, nelem);
    return 0;
  }

  (void) memcpy(array, value.data, value.size);

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
    (void) fprintf(stderr, "rtdb_seq_ma_get: handle (%d) is invalid\n", handle);
    return 0;
  }

  /* Retrieve the type info from the data base */

  if (!rtdb_seq_get_info(handle, name, ma_type, nelem, date)) {
    /*
    (void) fprintf(stderr, "rtdb_seq_ma_get: get info failed for \"%s\" in %s\n",
		   name, rtdb[handle].filename);
                   */
    return 0;
  }

  /* Allocate MA pointers */

  ma_type_buf = (Integer) *ma_type;
  nelem_buf   = (Integer) *nelem;

  if (!MA_allocate_heap(ma_type_buf, nelem_buf, name, &ma_handle_buf)) {
    (void) fprintf(stderr, "rtdb_seq_ma_get: MA_allocate_heap failed, nelem=%d\n",
		   *nelem);
    return 0;
  }
  *ma_handle = ma_handle_buf;
  
  if (!MA_get_pointer(ma_handle_buf, &ma_data)) {
    (void) fprintf(stderr, "rtdb_seq_ma_get: MA_get_pointer failed\n");
    return 0;
  }

  if (!rtdb_seq_get(handle, name, *ma_type, *nelem, ma_data)) {
    (void) fprintf(stderr, "rtdb_seq_ma_get: rtdb_seq_get failed for %s\n", name);
    (void) MA_free_heap(ma_handle_buf);
    return 0;
  }

  return 1;
}
