/*$Id$*/
#include <stdio.h>
#include <string.h>
#include "rtdb.h"
#include "macdecls.h"
#include "ga.h"
#include "hdbm/hdbm.h"

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

/* Mode in which DB was opended */
static int par_mode[MAX_RTDB] = {INACTIVE,INACTIVE,INACTIVE,INACTIVE,INACTIVE};
static int parallel_mode=RTDB_PAR_MODE; /* Current mode RTDB_SEQ_MODE/RTDB_PAR_MODE */

/**
\ingroup rtdb 
@{
*/

/**
  \brief Change the RTDB access mode

  The RTDB can be accessed in a serial or parallel mode. In either case only
  process rank 0 actually interacts with the RTDB. In parallel mode accesses
  to the RTDB behave as collective operations. I.e. all processes block on 
  rtdb_put operations in addition on rtdb_get operations all processes receive
  the same data. In serial mode only process rank 0 should make RTDB calls.

  \param mode [Input] the new access mode (valid values RTDB_SEQ_MODE, RTDB_PAR_MODE)

  \return the old access mode value
*/
int rtdb_parallel(const int mode)
/*
  Set the parallel access mode of all databases to mode and
  return the previous setting
*/
{
  int old = parallel_mode;

  if (mode)
    parallel_mode = RTDB_PAR_MODE;
  else
    parallel_mode = RTDB_SEQ_MODE;
  return old;
}

static int verify_parallel_access()
/*
  Return true if access mode / processor values are sensible
*/
{
#ifdef GAGROUPS
  int me = GA_Nodeid();
#endif
  if ((parallel_mode == RTDB_SEQ_MODE) && (me != 0)) {
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
  Integer len = MA_sizeof(ma_type, nelem, MT_CHAR);
  Integer from = 0;
  GA_Brdcst(data, len, from);
}

/**
  \brief Open an RTDB stored on a given file

  The RTDB is stored on a file. Before the RTDB can be accessed in the program
  the file needs to be opened. To access the RTDB a handle is associated with
  with the file and this handle is used in the actual access routines.

  \param filename [Input] the name of the file holding the RTDB
  \param mode     [Input] the initial access mode of the RTDB
  \param handle   [Output] the RTBD handle

  \return The return value of the file open command
*/
int rtdb_open(const char *filename, const char *mode, int *handle)
{
  int status;
#ifdef GAGROUPS
  int me = GA_Nodeid();
#else
  me = GA_Nodeid();
#endif

  if (!verify_parallel_access()) return 0;

  if (parallel_mode == RTDB_SEQ_MODE || me == 0)
    status = rtdb_seq_open(filename, mode, handle);

  if (parallel_mode == RTDB_PAR_MODE) {
    rtdb_broadcast(TYPE_RTDB_HANDLE, MT_INT, 1, (void *) handle);
    rtdb_broadcast(TYPE_RTDB_STATUS, MT_INT, 1, (void *) &status);
  }

  if (status)
    par_mode[*handle] = parallel_mode;

  return status;
}

/**
  \brief Clone the RTDB creating a new file

  Take the filename from the specified RTDB and create a new filename with
  the specified suffix. Then copy the current RTDB to the new file.

  \param handle [Input] the RTDB handle
  \param suffix [Input] the suffix for the new RTDB file

  \return the status of the file copy operation
*/
int rtdb_clone(const int handle, const char *suffix)
{
  int status;
#ifdef GAGROUPS
  int me = GA_Nodeid();
#else
  me = GA_Nodeid();
#endif

  if (!verify_parallel_access()) return 0;
  if (handle < 0 || handle >= MAX_RTDB) {
    (void) fprintf(stderr, "rtdb_clone: handle out of range %d\n", handle);
    (void) fflush(stderr);
    return 0;
  }
  if (par_mode[handle] == INACTIVE) {
    (void) fprintf(stderr, "rtdb_clone: handle not active %d\n",handle);
    (void) fflush(stderr);
    return 0;
  }
  if (parallel_mode != par_mode[handle]) {
    (void) fprintf(stderr, "rtdb_clone: mode of open and copy mismatch\n");
    (void) fflush(stderr);
    return 0;
  }

  if (parallel_mode == RTDB_SEQ_MODE || me == 0)
    status = rtdb_seq_copy(handle, suffix);

  if (parallel_mode == RTDB_PAR_MODE) 
    rtdb_broadcast(TYPE_RTDB_STATUS, MT_INT, 1, (void *) &status);

  return status;
}
int rtdb_getfname(const int handle,
		  char* fname)
{
  int status;
  int length;
#ifdef GAGROUPS
  int me = GA_Nodeid();
#endif

  if (!verify_parallel_access()) return 0;

  if (handle < 0 || handle >= MAX_RTDB) {
    (void) fprintf(stderr, "rtdb_getfname: handle out of range %d\n", handle);
    (void) fflush(stderr);
    return 0;
  }
  if (par_mode[handle] == INACTIVE) {
    (void) fprintf(stderr, "rtdb_getfname: handle not active %d\n",handle);
    (void) fflush(stderr);
    return 0;
  }

  if (par_mode[handle] == RTDB_SEQ_MODE && parallel_mode == RTDB_PAR_MODE) {
    (void) fprintf(stderr, "rtdb_getfname: seq. open and par. get\n");
    (void) fflush(stderr);
    return 0;
  }

  if (parallel_mode == RTDB_SEQ_MODE || me == 0)
    status = rtdb_seq_getfname(handle, fname);

  if (parallel_mode == RTDB_PAR_MODE) {
    if ( me == 0) length = strlen(fname)+1;
    rtdb_broadcast(TYPE_RTDB_STATUS, MT_INT, 1, (void *) &status);
    rtdb_broadcast(TYPE_RTDB_STATUS, MT_INT, 1, (void *) &length);
    
    if (status) {
      rtdb_broadcast(TYPE_RTDB_FNAME,  MT_CHAR, length,  (void *) fname);
    }
  }

  return status;
}

/**
  \brief Close the RTDB

  Close the specified RTDB, the handle is no longer valid after this operation.
  \param handle [Input] the RTDB handle
  \param mode   [Input] the access mode on close, this has to match the actual
  current access mode

  \return the status from the file close operation
*/
int rtdb_close(const int handle, const char *mode)
{
  int status;
#ifdef GAGROUPS
  int me = GA_Nodeid();
#endif
  
  if (!verify_parallel_access()) return 0;
  if (handle < 0 || handle >= MAX_RTDB) {
    (void) fprintf(stderr, "rtdb_close: handle out of range %d\n", handle);
    (void) fflush(stderr);
    return 0;
  }
  if (par_mode[handle] == INACTIVE) {
    (void) fprintf(stderr, "rtdb_close: handle not active %d\n",handle);
    (void) fflush(stderr);
    return 0;
  }
  
  if (parallel_mode != par_mode[handle]) {
    (void) fprintf(stderr, "rtdb_close: mode of open and close mismatch\n");
    (void) fflush(stderr);
    return 0;
  }

  if (par_mode[handle] == RTDB_SEQ_MODE || me == 0)
    status = rtdb_seq_close(handle, mode);

  if (parallel_mode == RTDB_PAR_MODE)
    rtdb_broadcast(TYPE_RTDB_STATUS, MT_INT, 1, (void *) &status);

  if (status) 
    par_mode[handle] = INACTIVE;

  return status;
}

/**
  \brief Store data on the RTDB

  Store data on the specified RTDB. The data is stored with a key by which it
  can be reference for retrieval.

  \param handle  [Input] the RTDB handle
  \param name    [Input] the key for the data
  \param ma_type [Input] the type of the data specified by one of the MA data types (see mafdecls.fh)
  \param nelem   [Input] the number of elements of the specified type
  \param array   [Input] the actual data

  \return the status of the write operation
*/
int rtdb_put(const int handle, const char *name, const int ma_type,
	     const int nelem, const void *array)
{
  int status;
#ifdef GAGROUPS
  int me = GA_Nodeid();
#endif

  if (!verify_parallel_access()) return 0;

  if (handle < 0 || handle >= MAX_RTDB) {
    (void) fprintf(stderr, "rtdb_put: handle out of range %d\n", handle);
    (void) fflush(stderr);
    return 0;
  }
  if (par_mode[handle] == INACTIVE) {
    (void) fprintf(stderr, "rtdb_put: handle not active %d\n",handle);
    (void) fflush(stderr);
    return 0;
  }

  if (par_mode[handle] == RTDB_SEQ_MODE && parallel_mode == RTDB_PAR_MODE) {
    (void) fprintf(stderr, "rtdb_put: seq. open and par. put\n");
    (void) fflush(stderr);
    return 0;
  }

  if (parallel_mode == RTDB_SEQ_MODE || me == 0)
    status = rtdb_seq_put(handle, name, ma_type, nelem, array);
    
  if (parallel_mode == RTDB_PAR_MODE)
    rtdb_broadcast(TYPE_RTDB_STATUS, MT_INT, 1, (void *) &status);


  return status;
}

/**
  \brief Retrieve data from the RTDB

  Retrieve the data associated with the specified key from the RTDB and 
  return it in the provided array.

  \param handle  [Input] the RTDB handle
  \param name    [Input] the key for the data
  \param ma_type [Input] the type of the data specified by one of the MA data types (see mafdecls.fh)
  \param nelem   [Input] the number of elements of the specified type
  \param array   [Output] the actual data retrieved

  \return the status of the read operation
*/
int rtdb_get(const int handle, const char *name, const int ma_type,
		 const int nelem, void *array)
{
  int status;
#ifdef GAGROUPS
  int me = GA_Nodeid();
#endif

  if (!verify_parallel_access()) return 0;

  if (handle < 0 || handle >= MAX_RTDB) {
    (void) fprintf(stderr, "rtdb_get: handle out of range %d\n", handle);
    (void) fflush(stderr);
    return 0;
  }
  if (par_mode[handle] == INACTIVE) {
    (void) fprintf(stderr, "rtdb_get: handle not active %d\n",handle);
    (void) fflush(stderr);
    return 0;
  }

  if (par_mode[handle] == RTDB_SEQ_MODE && parallel_mode == RTDB_PAR_MODE) {
    (void) fprintf(stderr, "rtdb_get: seq. open and par. get\n");
    (void) fflush(stderr);
    return 0;
  }

  if (parallel_mode == RTDB_SEQ_MODE || me == 0)
    status = rtdb_seq_get(handle, name, ma_type, nelem, array);

  /* Implicit assumption that all processes call this with nelem
     greater than or equal to value used on process 0 */

  if (parallel_mode == RTDB_PAR_MODE) {
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
#ifdef GAGROUPS
  int me = GA_Nodeid();
#endif

  if (!verify_parallel_access()) return 0;

  if (handle < 0 || handle >= MAX_RTDB) {
    (void) fprintf(stderr, "rtdb_get_info: handle out of range %d\n", handle);
    (void) fflush(stderr);
    return 0;
  }
  if (par_mode[handle] == INACTIVE) {
    (void) fprintf(stderr, "rtdb_get_info: handle not active %d\n",handle);
    (void) fflush(stderr);
    return 0;
  }

  if (par_mode[handle] == RTDB_SEQ_MODE && parallel_mode == RTDB_PAR_MODE) {
    (void) fprintf(stderr, "rtdb_get_info: seq. open and par. get\n");
    (void) fflush(stderr);
    return 0;
  }

  if (parallel_mode == RTDB_SEQ_MODE || me == 0)
    status = rtdb_seq_get_info(handle, name, ma_type, nelem, date);

  if (parallel_mode == RTDB_PAR_MODE) {
    rtdb_broadcast(TYPE_RTDB_STATUS, MT_INT, 1, (void *) &status);
    
    if (status) {
      rtdb_broadcast(TYPE_RTDB_NELEM, MT_INT,  1,   (void *) nelem);
      rtdb_broadcast(TYPE_RTDB_TYPE,  MT_INT,  1,   (void *) ma_type);
      rtdb_broadcast(TYPE_RTDB_DATE,  MT_CHAR, 26,  (void *) date);
    }
  }

  return status;
}

/**
  \brief Retrieve data from the RTDB of unknown size

  Retrieve the data associated with the specified key from the RTDB. 
  In this case the size of the data is not known a priory. Hence the size is
  retrieved from the RTDB as well, an array of the appropriate size is
  allocated, and the MA handle as well as the size are returned. The calling
  program is responsible for deallocating the memory.

  \param handle  [Input] the RTDB handle
  \param name    [Input] the key for the data
  \param ma_type [Output] the type of the data specified by one of the MA data types (see mafdecls.fh)
  \param nelem   [Output] the number of elements of the specified type
  \param array   [Output] the actual data retrieved

  \return the status of the read operation
*/
int rtdb_ma_get(const int handle, const char *name, int *ma_type,
		    int *nelem, int *ma_handle)
{
  int status;
#ifdef GAGROUPS
  int me = GA_Nodeid();
#endif

  if (!verify_parallel_access()) return 0;

  if (handle < 0 || handle >= MAX_RTDB) {
    (void) fprintf(stderr, "rtdb_ma_get: handle out of range %d\n", handle);
    (void) fflush(stderr);
    return 0;
  }
  if (par_mode[handle] == INACTIVE) {
    (void) fprintf(stderr, "rtdb_ma_get: handle not active %d\n",handle);
    (void) fflush(stderr);
    return 0;
  }

  if (par_mode[handle] == RTDB_SEQ_MODE && parallel_mode == RTDB_PAR_MODE) {
    (void) fprintf(stderr, "rtdb_ma_get: seq. open and par. get\n");
    (void) fflush(stderr);
    return 0;
  }

  if (parallel_mode == RTDB_SEQ_MODE || me == 0)
    status = rtdb_seq_ma_get(handle, name, ma_type, nelem, ma_handle);

  if (parallel_mode == RTDB_PAR_MODE) {
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

	  GA_Error("rtdb_get_ma: rtdb_get_ma: MA_alloc failed, nelem=", (Integer)nelem);
	}
        *ma_handle = (int) ma_handle_buf;
      }
      
      if (!MA_get_pointer((Integer) *ma_handle, &ma_data)) {
	GA_Error("rtdb_get_ma: rtdb_get_ma: MA_get_ptr failed, nelem=", (Integer) nelem);
      }

      rtdb_broadcast(TYPE_RTDB_ARRAY, *ma_type, *nelem, (void *) ma_data);
    }
  }

  return status;
}

int rtdb_first(const int handle, const int namelen, char *name)
{
  int status;
#ifdef GAGROUPS
  int me = GA_Nodeid();
#endif

  if (!verify_parallel_access()) return 0;

  if (handle < 0 || handle >= MAX_RTDB) {
    (void) fprintf(stderr, "rtdb_first: handle out of range %d\n", handle);
    (void) fflush(stderr);
    return 0;
  }
  if (par_mode[handle] == INACTIVE) {
    (void) fprintf(stderr, "rtdb_first: handle not active %d\n",handle);
    (void) fflush(stderr);
    return 0;
  }

  if (par_mode[handle] == RTDB_SEQ_MODE && parallel_mode == RTDB_PAR_MODE) {
    (void) fprintf(stderr, "rtdb_first: seq. open and par. first\n");
    (void) fflush(stderr);
    return 0;
  }

  if (parallel_mode == RTDB_SEQ_MODE || me == 0)
    status = rtdb_seq_first(handle, namelen, name);

  if (parallel_mode == RTDB_PAR_MODE) {
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
#ifdef GAGROUPS
  int me = GA_Nodeid();
#endif

  if (!verify_parallel_access()) return 0;

  if (handle < 0 || handle >= MAX_RTDB) {
    (void) fprintf(stderr, "rtdb_next: handle out of range %d\n", handle);
    (void) fflush(stderr);
    return 0;
  }
  if (par_mode[handle] == INACTIVE) {
    (void) fprintf(stderr, "rtdb_next: handle not active %d\n",handle);
    (void) fflush(stderr);
    return 0;
  }

  if (par_mode[handle] == RTDB_SEQ_MODE && parallel_mode == RTDB_PAR_MODE) {
    (void) fprintf(stderr, "rtdb_next: seq. open and par. next\n");
    (void) fflush(stderr);
    return 0;
  }

  if (parallel_mode == RTDB_SEQ_MODE || me == 0)
    status = rtdb_seq_next(handle, namelen, name);
  
  if (parallel_mode == RTDB_PAR_MODE) {
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
/**
  \brief Print contents of the RTDB

  Prints the contents of the specified RTDB. An additional flag specifies
  whether only the keys should be printed, or the keys and their values.

  \param handle [Input] the RTDB handle
  \param print_values [Input] if TRUE print the keys and the values, otherwise
  just print the keys

  \return the status of the print operation
*/
int rtdb_print(const int handle, const int print_values)
{
  int status;
#ifdef GAGROUPS
  int me = GA_Nodeid();
#endif

  if (!verify_parallel_access()) return 0;

  if (handle < 0 || handle >= MAX_RTDB) {
    (void) fprintf(stderr, "rtdb_print: handle out of range %d\n", handle);
    (void) fflush(stderr);
    return 0;
  }
  if (par_mode[handle] == INACTIVE) {
    (void) fprintf(stderr, "rtdb_print: handle not active %d\n",handle);
    (void) fflush(stderr);
    return 0;
  }

  if (par_mode[handle] == RTDB_SEQ_MODE && parallel_mode == RTDB_PAR_MODE) {
    (void) fprintf(stderr, "rtdb_print: seq. open and par. print\n");
    (void) fflush(stderr);
    return 0;
  }

  if (parallel_mode == RTDB_SEQ_MODE || me == 0)
    status = rtdb_seq_print(handle, print_values);

  if (parallel_mode == RTDB_PAR_MODE)
    rtdb_broadcast(TYPE_RTDB_STATUS, MT_INT, 1, (void *) &status);

  return status;
}

/**
  \brief Delete the data associated with a key from the RTDB

  \param handle [Input] the RTDB handle
  \param name   [Input] the key 

  \return the status of the delete operation
*/
int rtdb_delete(const int handle, const char *name)
{
  int status;
#ifdef GAGROUPS
  int me = GA_Nodeid();
#endif

  if (!verify_parallel_access()) return 0;

  if (handle < 0 || handle >= MAX_RTDB) {
    (void) fprintf(stderr, "rtdb_delete: handle out of range %d\n", handle);
    (void) fflush(stderr);
    return 0;
  }
  if (par_mode[handle] == INACTIVE) {
    (void) fprintf(stderr, "rtdb_delete: handle not active %d\n",handle);
    (void) fflush(stderr);
    return 0;
  }

  if (par_mode[handle] == RTDB_SEQ_MODE && parallel_mode == RTDB_PAR_MODE) {
    (void) fprintf(stderr, "rtdb_delete: seq. open and par. delete\n");
    (void) fflush(stderr);
    return 0;
  }

  if (parallel_mode == RTDB_SEQ_MODE || me == 0)
    status = rtdb_seq_delete(handle, name);

  if (parallel_mode == RTDB_PAR_MODE)
    rtdb_broadcast(TYPE_RTDB_STATUS, MT_INT, 1, (void *) &status);

  return status;
}

  
void rtdb_print_usage()
{
#ifdef USE_HDBM
#ifdef GAGROUPS
  int me = GA_Nodeid();
#endif
  if (me == 0)
    hdbm_print_usage();
#endif
}

/**
@}
*/
