/*$Id$*/
#include <stdio.h>
#include <string.h>
#include "srtdb.h"
#include "macdecls.h"
#include "global.h"
#include "hdbm.h"

#define TYPE_RTDB_HANDLE  30001
#define TYPE_RTDB_STATUS  30002
#define TYPE_RTDB_NELEM   30003 
#define TYPE_RTDB_ARRAY   30004 
#define TYPE_RTDB_LEN     30005
#define TYPE_RTDB_NAME    30006
#define TYPE_RTDB_DATE    30007
#define TYPE_RTDB_TYPE    30008
#define TYPE_RTDB_FNAME   30009
/*#define GAGROUPS 1*/

#ifndef GAGROUPS
#endif

#define MAX_RTDB 5
#define INACTIVE  -1
#define SEQUENTIAL 0
#define PARALLEL   1

int srtdb_open(const char *filename, const char *mode, int *handle)
{
  int status;

    status = srtdb_seq_open(filename, mode, handle);

  return status;
}



int srtdb_close(const int handle, const char *mode)
{
  int status;
  
  if (handle < 0 || handle >= MAX_RTDB) {
    (void) fprintf(stderr, "srtdb_close: handle out of range %d\n", handle);
    (void) fflush(stderr);
    return 0;
  }
  

    status = srtdb_seq_close(handle, mode);

  return status;
}

int srtdb_put(const int handle, const char *name, const int ma_type,
	     const int nelem, const void *array)
{
  int status;

  if (handle < 0 || handle >= MAX_RTDB) {
    (void) fprintf(stderr, "srtdb_put: handle out of range %d\n", handle);
    (void) fflush(stderr);
    return 0;
  }


    status = srtdb_seq_put(handle, name, ma_type, nelem, array);
    

  return status;
}

int srtdb_get(const int handle, const char *name, const int ma_type,
		 const int nelem, void *array)
{
  int status;

  if (handle < 0 || handle >= MAX_RTDB) {
    (void) fprintf(stderr, "srtdb_get: handle out of range %d\n", handle);
    (void) fflush(stderr);
    return 0;
  }

    status = srtdb_seq_get(handle, name, ma_type, nelem, array);

  return status;
}

int srtdb_get_info(const int handle,
		  const char *name, int *ma_type, int *nelem, char date[26])
{
  int status;

  if (handle < 0 || handle >= MAX_RTDB) {
    (void) fprintf(stderr, "srtdb_get_info: handle out of range %d\n", handle);
    (void) fflush(stderr);
    return 0;
  }

    status = srtdb_seq_get_info(handle, name, ma_type, nelem, date);

  return status;
}

int srtdb_first(const int handle, const int namelen, char *name)
{
  int status;

  if (handle < 0 || handle >= MAX_RTDB) {
    (void) fprintf(stderr, "srtdb_first: handle out of range %d\n", handle);
    (void) fflush(stderr);
    return 0;
  }


    status = srtdb_seq_first(handle, namelen, name);

  return status;
}

int srtdb_next(const int handle, const int namelen, char *name)
{
  int status;

  if (handle < 0 || handle >= MAX_RTDB) {
    (void) fprintf(stderr, "srtdb_next: handle out of range %d\n", handle);
    (void) fflush(stderr);
    return 0;
  }

    status = srtdb_seq_next(handle, namelen, name);
  
  return status;
}

int srtdb_print(const int handle, const int print_values)
{
  int status;

  if (handle < 0 || handle >= MAX_RTDB) {
    (void) fprintf(stderr, "srtdb_print: handle out of range %d\n", handle);
    (void) fflush(stderr);
    return 0;
  }

    status = srtdb_seq_print(handle, print_values);

  return status;
}

int srtdb_delete(const int handle, const char *name)
{
  int status;

  if (handle < 0 || handle >= MAX_RTDB) {
    (void) fprintf(stderr, "srtdb_delete: handle out of range %d\n", handle);
    (void) fflush(stderr);
    return 0;
  }

    status = srtdb_seq_delete(handle, name);

  return status;
}

  
void srtdb_print_usage()
{
   hdbm_print_usage();
}
