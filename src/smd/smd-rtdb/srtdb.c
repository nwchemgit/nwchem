/*$Id: srtdb.c,v 1.1 2008-04-21 18:19:27 marat Exp $*/
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
static int me;
#endif

#define MAX_RTDB 5
#define INACTIVE  -1
#define SEQUENTIAL 0
#define PARALLEL   1

/* Mode in which DB was opended */
static int par_mode[MAX_RTDB] = {INACTIVE,INACTIVE,INACTIVE,INACTIVE,INACTIVE};
static int parallel_mode=PARALLEL; /* Current mode SEQUENTIAL/PARALLEL */

int srtdb_parallel(const int mode)
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
#ifdef GAGROUPS
  int me = ga_nodeid_();
#endif
  if ((parallel_mode == SEQUENTIAL) && (me != 0)) {
    (void) fflush(stdout);
    (void) fprintf(stderr,"rtdb: sequential access only possible for process 0\n");
    (void) fflush(stderr);
    return 0;
  }
  else
    return 1;
}

static void srtdb_broadcast(const int msg_type, const int ma_type, 
			   const int nelem, void *data)
/*
  Convenience routine for rtdb that broadcasts MA data types from
  node 0 to all other nodes
*/
{
  Integer len = MA_sizeof(ma_type, nelem, MT_CHAR);
  Integer from = 0;
  Integer type = msg_type;
  ga_brdcst_(&type, (char *) data, &len, &from);
}

int srtdb_open(const char *filename, const char *mode, int *handle)
{
  int status;
#ifdef GAGROUPS
  int me = ga_nodeid_();
#else
  me = ga_nodeid_();
#endif

  if (!verify_parallel_access()) return 0;

  if (parallel_mode == SEQUENTIAL || me == 0)
    status = srtdb_seq_open(filename, mode, handle);

  if (parallel_mode == PARALLEL) {
    srtdb_broadcast(TYPE_RTDB_HANDLE, MT_INT, 1, (void *) handle);
    srtdb_broadcast(TYPE_RTDB_STATUS, MT_INT, 1, (void *) &status);
  }

  if (status)
    par_mode[*handle] = parallel_mode;

  return status;
}

int srtdb_clone(const int handle, const char *suffix)
{
  int status;
#ifdef GAGROUPS
  int me = ga_nodeid_();
#else
  me = ga_nodeid_();
#endif

  if (!verify_parallel_access()) return 0;
  if (handle < 0 || handle >= MAX_RTDB) {
    (void) fprintf(stderr, "srtdb_clone: handle out of range %d\n", handle);
    (void) fflush(stderr);
    return 0;
  }
  if (parallel_mode != par_mode[handle]) {
    (void) fprintf(stderr, "srtdb_clone: mode of open and copy mismatch\n");
    (void) fflush(stderr);
    return 0;
  }

  if (parallel_mode == SEQUENTIAL || me == 0)
    status = srtdb_seq_copy(handle, suffix);

  if (parallel_mode == PARALLEL) 
    srtdb_broadcast(TYPE_RTDB_STATUS, MT_INT, 1, (void *) &status);

  return status;
}
int srtdb_getfname(const int handle,
		  char fname[36])
{
  int status;
#ifdef GAGROUPS
  int me = ga_nodeid_();
#endif

  if (!verify_parallel_access()) return 0;

  if (handle < 0 || handle >= MAX_RTDB) {
    (void) fprintf(stderr, "srtdb_getfname: handle out of range %d\n", handle);
    (void) fflush(stderr);
    return 0;
  }

  if (par_mode[handle] == SEQUENTIAL && parallel_mode == PARALLEL) {
    (void) fprintf(stderr, "srtdb_getfname: seq. open and par. get\n");
    (void) fflush(stderr);
    return 0;
  }

  if (parallel_mode == SEQUENTIAL || me == 0)
    status = srtdb_seq_getfname(handle, fname);

  if (parallel_mode == PARALLEL) {
    srtdb_broadcast(TYPE_RTDB_STATUS, MT_INT, 1, (void *) &status);
    
    if (status) {
      srtdb_broadcast(TYPE_RTDB_FNAME,  MT_CHAR, 36,  (void *) fname);
    }
  }

  return status;
}

int srtdb_close(const int handle, const char *mode)
{
  int status;
#ifdef GAGROUPS
  int me = ga_nodeid_();
#endif
  
  if (!verify_parallel_access()) return 0;
  if (handle < 0 || handle >= MAX_RTDB) {
    (void) fprintf(stderr, "srtdb_close: handle out of range %d\n", handle);
    (void) fflush(stderr);
    return 0;
  }
  
  if (parallel_mode != par_mode[handle]) {
    (void) fprintf(stderr, "srtdb_close: mode of open and close mismatch\n");
    (void) fflush(stderr);
    return 0;
  }

  if (par_mode[handle] == SEQUENTIAL || me == 0)
    status = srtdb_seq_close(handle, mode);

  if (parallel_mode == PARALLEL)
    srtdb_broadcast(TYPE_RTDB_STATUS, MT_INT, 1, (void *) &status);
  return status;
}

int srtdb_put(const int handle, const char *name, const int ma_type,
	     const int nelem, const void *array)
{
  int status;
#ifdef GAGROUPS
  int me = ga_nodeid_();
#endif

  if (!verify_parallel_access()) return 0;

  if (handle < 0 || handle >= MAX_RTDB) {
    (void) fprintf(stderr, "srtdb_put: handle out of range %d\n", handle);
    (void) fflush(stderr);
    return 0;
  }

  if (par_mode[handle] == SEQUENTIAL && parallel_mode == PARALLEL) {
    (void) fprintf(stderr, "srtdb_put: seq. open and par. put\n");
    (void) fflush(stderr);
    return 0;
  }

  if (parallel_mode == SEQUENTIAL || me == 0)
    status = srtdb_seq_put(handle, name, ma_type, nelem, array);
    
  if (parallel_mode == PARALLEL)
    srtdb_broadcast(TYPE_RTDB_STATUS, MT_INT, 1, (void *) &status);


  return status;
}

int srtdb_get(const int handle, const char *name, const int ma_type,
		 const int nelem, void *array)
{
  int status;
#ifdef GAGROUPS
  int me = ga_nodeid_();
#endif

  if (!verify_parallel_access()) return 0;

  if (handle < 0 || handle >= MAX_RTDB) {
    (void) fprintf(stderr, "srtdb_get: handle out of range %d\n", handle);
    (void) fflush(stderr);
    return 0;
  }

  if (par_mode[handle] == SEQUENTIAL && parallel_mode == PARALLEL) {
    (void) fprintf(stderr, "srtdb_get: seq. open and par. get\n");
    (void) fflush(stderr);
    return 0;
  }

  if (parallel_mode == SEQUENTIAL || me == 0)
    status = srtdb_seq_get(handle, name, ma_type, nelem, array);

  /* Implicit assumption that all processes call this with nelem
     greater than or equal to value used on process 0 */

  if (parallel_mode == PARALLEL) {
    srtdb_broadcast(TYPE_RTDB_STATUS, MT_INT, 1, (void *) &status);
    
    if (status) {
      srtdb_broadcast(TYPE_RTDB_NELEM, MT_INT, 1, (void *) &nelem);
      srtdb_broadcast(TYPE_RTDB_ARRAY, ma_type, nelem, (void *) array);
    }
  }

  return status;
}

int srtdb_get_info(const int handle,
		  const char *name, int *ma_type, int *nelem, char date[26])
{
  int status;
#ifdef GAGROUPS
  int me = ga_nodeid_();
#endif

  if (!verify_parallel_access()) return 0;

  if (handle < 0 || handle >= MAX_RTDB) {
    (void) fprintf(stderr, "srtdb_get_info: handle out of range %d\n", handle);
    (void) fflush(stderr);
    return 0;
  }

  if (par_mode[handle] == SEQUENTIAL && parallel_mode == PARALLEL) {
    (void) fprintf(stderr, "srtdb_get_info: seq. open and par. get\n");
    (void) fflush(stderr);
    return 0;
  }

  if (parallel_mode == SEQUENTIAL || me == 0)
    status = srtdb_seq_get_info(handle, name, ma_type, nelem, date);

  if (parallel_mode == PARALLEL) {
    srtdb_broadcast(TYPE_RTDB_STATUS, MT_INT, 1, (void *) &status);
    
    if (status) {
      srtdb_broadcast(TYPE_RTDB_NELEM, MT_INT,  1,   (void *) nelem);
      srtdb_broadcast(TYPE_RTDB_TYPE,  MT_INT,  1,   (void *) ma_type);
      srtdb_broadcast(TYPE_RTDB_DATE,  MT_CHAR, 26,  (void *) date);
    }
  }

  return status;
}

int srtdb_ma_get(const int handle, const char *name, int *ma_type,
		    int *nelem, int *ma_handle)
{
  int status;
#ifdef GAGROUPS
  int me = ga_nodeid_();
#endif

  if (!verify_parallel_access()) return 0;

  if (handle < 0 || handle >= MAX_RTDB) {
    (void) fprintf(stderr, "srtdb_ma_get: handle out of range %d\n", handle);
    (void) fflush(stderr);
    return 0;
  }

  if (par_mode[handle] == SEQUENTIAL && parallel_mode == PARALLEL) {
    (void) fprintf(stderr, "srtdb_ma_get: seq. open and par. get\n");
    (void) fflush(stderr);
    return 0;
  }

  if (parallel_mode == SEQUENTIAL || me == 0)
    status = srtdb_seq_ma_get(handle, name, ma_type, nelem, ma_handle);

  if (parallel_mode == PARALLEL) {
    srtdb_broadcast(TYPE_RTDB_STATUS, MT_INT, 1, (void *) &status);
    
    if (status) {
      void *ma_data;
      
      srtdb_broadcast(TYPE_RTDB_TYPE,  MT_INT, 1, (void *) ma_type);
      srtdb_broadcast(TYPE_RTDB_NELEM, MT_INT, 1, (void *) nelem);

      if (me > 0) {
        Integer ma_handle_buf;
        Integer ma_type_buf = *ma_type;
        Integer nelem_buf = *nelem;
	if (!MA_allocate_heap(ma_type_buf, nelem_buf, name, &ma_handle_buf)) {
	  /* Cannot just return and error here since cannot propogate the
	     error condition to all the other processes ... exit via
	     the TCGMSG error routine */

	  ga_error("srtdb_get_ma: srtdb_get_ma: MA_alloc failed, nelem=",
		(Integer) nelem);
	}
        *ma_handle = (int) ma_handle_buf;
      }
      
      if (!MA_get_pointer((Integer) *ma_handle, &ma_data)) {
	ga_error("srtdb_get_ma: srtdb_get_ma: MA_get_ptr failed, nelem=",
	      (Integer) nelem);
      }

      srtdb_broadcast(TYPE_RTDB_ARRAY, *ma_type, *nelem, (void *) ma_data);
    }
  }

  return status;
}

int srtdb_first(const int handle, const int namelen, char *name)
{
  int status;
#ifdef GAGROUPS
  int me = ga_nodeid_();
#endif

  if (!verify_parallel_access()) return 0;

  if (handle < 0 || handle >= MAX_RTDB) {
    (void) fprintf(stderr, "srtdb_first: handle out of range %d\n", handle);
    (void) fflush(stderr);
    return 0;
  }

  if (par_mode[handle] == SEQUENTIAL && parallel_mode == PARALLEL) {
    (void) fprintf(stderr, "srtdb_first: seq. open and par. first\n");
    (void) fflush(stderr);
    return 0;
  }

  if (parallel_mode == SEQUENTIAL || me == 0)
    status = srtdb_seq_first(handle, namelen, name);

  if (parallel_mode == PARALLEL) {
    srtdb_broadcast(TYPE_RTDB_STATUS, MT_INT, 1, (void *) &status);
    
    if (status) {
      int len;
      if (me == 0)len = strlen(name)+1;
      srtdb_broadcast(TYPE_RTDB_LEN, MT_INT, 1, (void *) &len);
      srtdb_broadcast(TYPE_RTDB_NAME, MT_CHAR, len, (void *) name);
    }
  }
  
  return status;
}

int srtdb_next(const int handle, const int namelen, char *name)
{
  int status;
#ifdef GAGROUPS
  int me = ga_nodeid_();
#endif

  if (!verify_parallel_access()) return 0;

  if (handle < 0 || handle >= MAX_RTDB) {
    (void) fprintf(stderr, "srtdb_next: handle out of range %d\n", handle);
    (void) fflush(stderr);
    return 0;
  }

  if (par_mode[handle] == SEQUENTIAL && parallel_mode == PARALLEL) {
    (void) fprintf(stderr, "srtdb_next: seq. open and par. next\n");
    (void) fflush(stderr);
    return 0;
  }

  if (parallel_mode == SEQUENTIAL || me == 0)
    status = srtdb_seq_next(handle, namelen, name);
  
  if (parallel_mode == PARALLEL) {
    srtdb_broadcast(TYPE_RTDB_STATUS, MT_INT, 1, (void *) &status);
    
    if (status) {
      int len;
      if (me == 0)len = strlen(name)+1;
      srtdb_broadcast(TYPE_RTDB_LEN, MT_INT, 1, (void *) &len);
      srtdb_broadcast(TYPE_RTDB_NAME, MT_CHAR, len, (void *) name);
    }
    
  }

  return status;
}

int srtdb_print(const int handle, const int print_values)
{
  int status;
#ifdef GAGROUPS
  int me = ga_nodeid_();
#endif

  if (!verify_parallel_access()) return 0;

  if (handle < 0 || handle >= MAX_RTDB) {
    (void) fprintf(stderr, "srtdb_print: handle out of range %d\n", handle);
    (void) fflush(stderr);
    return 0;
  }

  if (par_mode[handle] == SEQUENTIAL && parallel_mode == PARALLEL) {
    (void) fprintf(stderr, "srtdb_print: seq. open and par. print\n");
    (void) fflush(stderr);
    return 0;
  }

  if (parallel_mode == SEQUENTIAL || me == 0)
    status = srtdb_seq_print(handle, print_values);

  if (parallel_mode == PARALLEL)
    srtdb_broadcast(TYPE_RTDB_STATUS, MT_INT, 1, (void *) &status);

  return status;
}

int srtdb_delete(const int handle, const char *name)
{
  int status;
#ifdef GAGROUPS
  int me = ga_nodeid_();
#endif

  if (!verify_parallel_access()) return 0;

  if (handle < 0 || handle >= MAX_RTDB) {
    (void) fprintf(stderr, "srtdb_delete: handle out of range %d\n", handle);
    (void) fflush(stderr);
    return 0;
  }

  if (par_mode[handle] == SEQUENTIAL && parallel_mode == PARALLEL) {
    (void) fprintf(stderr, "srtdb_delete: seq. open and par. delete\n");
    (void) fflush(stderr);
    return 0;
  }

  if (parallel_mode == SEQUENTIAL || me == 0)
    status = srtdb_seq_delete(handle, name);

  if (parallel_mode == PARALLEL)
    srtdb_broadcast(TYPE_RTDB_STATUS, MT_INT, 1, (void *) &status);

  return status;
}

  
void srtdb_print_usage()
{
#ifdef USE_HDBM
#ifdef GAGROUPS
  int me = ga_nodeid_();
#endif
  if (me == 0)
    hdbm_print_usage();
#endif
}
