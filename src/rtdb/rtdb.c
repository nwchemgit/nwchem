/*$Id: rtdb.c,v 1.13 1999-11-13 03:02:22 bjohnson Exp $*/
#include <stdio.h>
#include <string.h>
#include "rtdb.h"
#include "macdecls.h"
#include "sndrcv.h"
#include "global.h"

typedef long integer;		/* Equivalent C type to FORTRAN integer */

#define TYPE_RTDB_HANDLE  30001
#define TYPE_RTDB_STATUS  30002
#define TYPE_RTDB_NELEM   30003 
#define TYPE_RTDB_ARRAY   30004 
#define TYPE_RTDB_LEN     30005
#define TYPE_RTDB_NAME    30006
#define TYPE_RTDB_DATE    30007
#define TYPE_RTDB_TYPE    30008

static int me;

#define MAX_RTDB 5
#define INACTIVE  -1
#define SEQUENTIAL 0
#define PARALLEL   1

/* Mode in which DB was opended */
static int par_mode[MAX_RTDB] = {INACTIVE,INACTIVE,INACTIVE,INACTIVE,INACTIVE};
static int parallel_mode=PARALLEL; /* Current mode SEQUENTIAL/PARALLEL */

int rtdb_parallel(const int mode)
/*
  Set the parallel access mode of all databases to mode and
  return the previous setting
*/
{
  int old = parallel_mode;

  if (mode)
    parallel_mode = PARALLEL;
  else
    parallel_mode = SEQUENTIAL;
  
  return old;
}

static int verify_parallel_access()
/*
  Return true if access mode / processor values are sensible
*/
{
  if ((parallel_mode == SEQUENTIAL) && (me != 0)) {
    (void) fflush(stdout);
    (void) fprintf(stderr,"rtdb: sequential access only possible for process 0\n");
    (void) fflush(stderr);
    return 0;
  }
  else
    return 1;
}

static void rtdb_broadcast(const int msg_type, const int ma_type, 
			   const int nelem, void *data)
/*
  Convenience routine for rtdb that broadcasts MA data types from
  node 0 to all other nodes
*/
{
  long len = MA_sizeof(ma_type, nelem, MT_CHAR);
  long from = 0;
  long type = msg_type;
  ga_brdcst_(&type, (char *) data, &len, &from);
}

int rtdb_open(const char *filename, const char *mode, int *handle)
{
  int status;
  me = NODEID_();

  if (!verify_parallel_access()) return 0;

  if (parallel_mode == SEQUENTIAL || me == 0)
    status = rtdb_seq_open(filename, mode, handle);

  if (parallel_mode == PARALLEL) {
    rtdb_broadcast(TYPE_RTDB_HANDLE, MT_INT, 1, (void *) handle);
    rtdb_broadcast(TYPE_RTDB_STATUS, MT_INT, 1, (void *) &status);
  }

  if (status)
    par_mode[*handle] = parallel_mode;

  return status;
}

int rtdb_close(const int handle, const char *mode)
{
  int status;
  
  if (!verify_parallel_access()) return 0;
  if (handle < 0 || handle >= MAX_RTDB) {
    (void) fprintf(stderr, "rtdb_close: handle out of range %d\n", handle);
    (void) fflush(stderr);
    return 0;
  }
  
  if (parallel_mode != par_mode[handle]) {
    (void) fprintf(stderr, "rtdb_close: mode of open and close mismatch\n");
    (void) fflush(stderr);
    return 0;
  }

  if (par_mode[handle] == SEQUENTIAL || me == 0)
    status = rtdb_seq_close(handle, mode);

  if (parallel_mode == PARALLEL)
    rtdb_broadcast(TYPE_RTDB_STATUS, MT_INT, 1, (void *) &status);
  return status;
}

int rtdb_put(const int handle, const char *name, const int ma_type,
	     const int nelem, const void *array)
{
  int status;

  if (!verify_parallel_access()) return 0;

  if (handle < 0 || handle >= MAX_RTDB) {
    (void) fprintf(stderr, "rtdb_put: handle out of range %d\n", handle);
    (void) fflush(stderr);
    return 0;
  }

  if (par_mode[handle] == SEQUENTIAL && parallel_mode == PARALLEL) {
    (void) fprintf(stderr, "rtdb_put: seq. open and par. put\n");
    (void) fflush(stderr);
    return 0;
  }

  if (parallel_mode == SEQUENTIAL || me == 0)
    status = rtdb_seq_put(handle, name, ma_type, nelem, array);
    
  if (parallel_mode == PARALLEL)
    rtdb_broadcast(TYPE_RTDB_STATUS, MT_INT, 1, (void *) &status);


  return status;
}

int rtdb_get(const int handle, const char *name, const int ma_type,
		 const int nelem, void *array)
{
  int status;

  if (!verify_parallel_access()) return 0;

  if (handle < 0 || handle >= MAX_RTDB) {
    (void) fprintf(stderr, "rtdb_get: handle out of range %d\n", handle);
    (void) fflush(stderr);
    return 0;
  }

  if (par_mode[handle] == SEQUENTIAL && parallel_mode == PARALLEL) {
    (void) fprintf(stderr, "rtdb_get: seq. open and par. get\n");
    (void) fflush(stderr);
    return 0;
  }

  if (parallel_mode == SEQUENTIAL || me == 0)
    status = rtdb_seq_get(handle, name, ma_type, nelem, array);

  /* Implicit assumption that all processes call this with nelem
     greater than or equal to value used on process 0 */

  if (parallel_mode == PARALLEL) {
    rtdb_broadcast(TYPE_RTDB_STATUS, MT_INT, 1, (void *) &status);
    
    if (status) {
      rtdb_broadcast(TYPE_RTDB_NELEM, MT_INT, 1, (void *) &nelem);
      rtdb_broadcast(TYPE_RTDB_ARRAY, ma_type, nelem, (void *) array);
    }
  }

  return status;
}

int rtdb_get_info(const int handle,
		  const char *name, int *ma_type, int *nelem, char date[26])
{
  int status;

  if (!verify_parallel_access()) return 0;

  if (handle < 0 || handle >= MAX_RTDB) {
    (void) fprintf(stderr, "rtdb_get_info: handle out of range %d\n", handle);
    (void) fflush(stderr);
    return 0;
  }

  if (par_mode[handle] == SEQUENTIAL && parallel_mode == PARALLEL) {
    (void) fprintf(stderr, "rtdb_get_info: seq. open and par. get\n");
    (void) fflush(stderr);
    return 0;
  }

  if (parallel_mode == SEQUENTIAL || me == 0)
    status = rtdb_seq_get_info(handle, name, ma_type, nelem, date);

  if (parallel_mode == PARALLEL) {
    rtdb_broadcast(TYPE_RTDB_STATUS, MT_INT, 1, (void *) &status);
    
    if (status) {
      rtdb_broadcast(TYPE_RTDB_NELEM, MT_INT,  1,   (void *) nelem);
      rtdb_broadcast(TYPE_RTDB_TYPE,  MT_INT,  1,   (void *) ma_type);
      rtdb_broadcast(TYPE_RTDB_DATE,  MT_CHAR, 26,  (void *) date);
    }
  }

  return status;
}

int rtdb_ma_get(const int handle, const char *name, int *ma_type,
		    int *nelem, int *ma_handle)
{
  int status;

  if (!verify_parallel_access()) return 0;

  if (handle < 0 || handle >= MAX_RTDB) {
    (void) fprintf(stderr, "rtdb_ma_get: handle out of range %d\n", handle);
    (void) fflush(stderr);
    return 0;
  }

  if (par_mode[handle] == SEQUENTIAL && parallel_mode == PARALLEL) {
    (void) fprintf(stderr, "rtdb_ma_get: seq. open and par. get\n");
    (void) fflush(stderr);
    return 0;
  }

  if (parallel_mode == SEQUENTIAL || me == 0)
    status = rtdb_seq_ma_get(handle, name, ma_type, nelem, ma_handle);

  if (parallel_mode == PARALLEL) {
    rtdb_broadcast(TYPE_RTDB_STATUS, MT_INT, 1, (void *) &status);
    
    if (status) {
      void *ma_data;
      
      rtdb_broadcast(TYPE_RTDB_TYPE,  MT_INT, 1, (void *) ma_type);
      rtdb_broadcast(TYPE_RTDB_NELEM, MT_INT, 1, (void *) nelem);

      if (me > 0) {
        Integer ma_handle_buf;
        Integer ma_type_buf = *ma_type;
        Integer nelem_buf = *nelem;
	if (!MA_allocate_heap(ma_type_buf, nelem_buf, name, &ma_handle_buf)) {
	  /* Cannot just return and error here since cannot propogate the
	     error condition to all the other processes ... exit via
	     the TCGMSG error routine */

	  Error("rtdb_get_ma: rtdb_get_ma: MA_alloc failed, nelem=",
		(long) nelem);
	}
        *ma_handle = (int) ma_handle_buf;
      }
      
      if (!MA_get_pointer((Integer) *ma_handle, &ma_data)) {
	Error("rtdb_get_ma: rtdb_get_ma: MA_get_ptr failed, nelem=",
	      (long) nelem);
      }

      rtdb_broadcast(TYPE_RTDB_ARRAY, *ma_type, *nelem, (void *) ma_data);
    }
  }

  return status;
}

int rtdb_first(const int handle, const int namelen, char *name)
{
  int status;

  if (!verify_parallel_access()) return 0;

  if (handle < 0 || handle >= MAX_RTDB) {
    (void) fprintf(stderr, "rtdb_first: handle out of range %d\n", handle);
    (void) fflush(stderr);
    return 0;
  }

  if (par_mode[handle] == SEQUENTIAL && parallel_mode == PARALLEL) {
    (void) fprintf(stderr, "rtdb_first: seq. open and par. first\n");
    (void) fflush(stderr);
    return 0;
  }

  if (parallel_mode == SEQUENTIAL || me == 0)
    status = rtdb_seq_first(handle, namelen, name);

  if (parallel_mode == PARALLEL) {
    rtdb_broadcast(TYPE_RTDB_STATUS, MT_INT, 1, (void *) &status);
    
    if (status) {
      int len;
      if (me == 0)len = strlen(name)+1;
      rtdb_broadcast(TYPE_RTDB_LEN, MT_INT, 1, (void *) &len);
      rtdb_broadcast(TYPE_RTDB_NAME, MT_CHAR, len, (void *) name);
    }
  }
  
  return status;
}

int rtdb_next(const int handle, const int namelen, char *name)
{
  int status;

  if (!verify_parallel_access()) return 0;

  if (handle < 0 || handle >= MAX_RTDB) {
    (void) fprintf(stderr, "rtdb_next: handle out of range %d\n", handle);
    (void) fflush(stderr);
    return 0;
  }

  if (par_mode[handle] == SEQUENTIAL && parallel_mode == PARALLEL) {
    (void) fprintf(stderr, "rtdb_next: seq. open and par. next\n");
    (void) fflush(stderr);
    return 0;
  }

  if (parallel_mode == SEQUENTIAL || me == 0)
    status = rtdb_seq_next(handle, namelen, name);
  
  if (parallel_mode == PARALLEL) {
    rtdb_broadcast(TYPE_RTDB_STATUS, MT_INT, 1, (void *) &status);
    
    if (status) {
      int len;
      if (me == 0)len = strlen(name)+1;
      rtdb_broadcast(TYPE_RTDB_LEN, MT_INT, 1, (void *) &len);
      rtdb_broadcast(TYPE_RTDB_NAME, MT_CHAR, len, (void *) name);
    }
    
  }

  return status;
}

int rtdb_print(const int handle, const int print_values)
{
  int status;

  if (!verify_parallel_access()) return 0;

  if (handle < 0 || handle >= MAX_RTDB) {
    (void) fprintf(stderr, "rtdb_print: handle out of range %d\n", handle);
    (void) fflush(stderr);
    return 0;
  }

  if (par_mode[handle] == SEQUENTIAL && parallel_mode == PARALLEL) {
    (void) fprintf(stderr, "rtdb_print: seq. open and par. print\n");
    (void) fflush(stderr);
    return 0;
  }

  if (parallel_mode == SEQUENTIAL || me == 0)
    status = rtdb_seq_print(handle, print_values);

  if (parallel_mode == PARALLEL)
    rtdb_broadcast(TYPE_RTDB_STATUS, MT_INT, 1, (void *) &status);

  return status;
}

int rtdb_delete(const int handle, const char *name)
{
  int status;

  if (!verify_parallel_access()) return 0;

  if (handle < 0 || handle >= MAX_RTDB) {
    (void) fprintf(stderr, "rtdb_delete: handle out of range %d\n", handle);
    (void) fflush(stderr);
    return 0;
  }

  if (par_mode[handle] == SEQUENTIAL && parallel_mode == PARALLEL) {
    (void) fprintf(stderr, "rtdb_delete: seq. open and par. delete\n");
    (void) fflush(stderr);
    return 0;
  }

  if (parallel_mode == SEQUENTIAL || me == 0)
    status = rtdb_seq_delete(handle, name);

  if (parallel_mode == PARALLEL)
    rtdb_broadcast(TYPE_RTDB_STATUS, MT_INT, 1, (void *) &status);

  return status;
}

  
void rtdb_print_usage()
{
#ifdef USE_HDBM
  if (me == 0)
    hdbm_print_usage();
#endif
}
