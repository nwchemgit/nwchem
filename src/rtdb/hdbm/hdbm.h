/*
 $Id$
 */

#ifndef HASHTABLE_H
#define HASHTABLE_H

#ifndef TRUE
#define TRUE 1
#endif

#ifndef FALSE
#define FALSE 0
#endif

/*
  Routines that support a simple in core hash table.

  *SUMMARY*

  typedef struct datum {
    void *dptr;
    int dsize;
  } datum;

  typedef hdbm int;

  DATUMS_MATCH(a, b)

  *RECOMMENDATIONS* 

         Ensure that datums extracted from the table are always
         freed using datum_free().  This routine complements
         datum_make(), datum_duplicate() and datum_wrap().
*/

typedef struct datum {		/* The fundamental data type */
  void *dptr;			/* Ndbm uses char * */
  int dsize;
} datum;

typedef int hdbm;		/* Handle to refer to hashtable */

#define DATUMS_MATCH(a, b) (a.dsize == b.dsize && \
			   memcmp(a.dptr, b.dptr, (size_t) a.dsize) == 0)

/*
  Determine if db is a valid handle for an open hashtable.
*/
extern int hdbm_check(hdbm);

/*
  Make a datum that refers to the region pointed to
  by the arguments .. NO MEMORY IS ALLOCATED and unless
  the provided pointer was returned from malloc() this
  datum cannot be freed using datum_free().
*/
extern void datum_wrap(const void *, int, datum *);

/*
  Make a datum by duplicating the data described
  by the arguments. 

  Return  1 if successful
          0 if failed to allocate necessary memory
*/
extern int datum_make(const void *, int, datum *);

/*
  Duplicate the datum d, allocating space as required.  Return 1
  on success, 0 on failure.
*/
extern int datum_duplicate(datum, datum *);

/*
  Print summary of datum to stdout in hex
*/
extern void datum_print_hex(datum);

/*
  Prints datum to stdout assuming a null terminated string.
*/
extern void datum_print_string(datum);

/*
  Free the memory pointed to by d.dptr ... must have come from malloc.
*/
extern void datum_free(datum);

/*
  Open a database with given name. 

  If (use_file)
     Store all keys in memory.  Store all keys and values on disk.
     Name is interpreted as the path to a file.
  else
     Store all keys and values in memory.  Nothing on disk.
     Name serves merely to identify the database.

  On sucess return true and the value of the handle.

  On failure return false and an invalid handle.
*/
extern int hdbm_open(const char *, int, hdbm *);

/*
  Close this db ... eventually want to flush all entries to disk
  etc., but right now is just core resident so just call hdbm_clear
  and mark as inactive.

  Return 1 on sucess, 0 on failure.
*/
extern int hdbm_close(hdbm);


/*
  Delete any original content of fileto and then copy the 
  database in filefrom into fileto.  Both files are closed
  at the end of the operation.
*/
extern int hdbm_file_copy(const char *, const char *);

/*
  Flush buffered I/O to disk
*/
extern int hdbm_file_flush(hdbm);

/*
  Compress out dead space from database on disk which
  must be closed 
*/
extern int hdbm_file_compress(const char *);

/*
  Insert a new entry into the table without looking
  for duplicates.  The datums (key and value) are duplicated.

  Return 1 on sucess, 0 on failure.

*/
extern int hdbm_insert(hdbm, datum, datum);

/*
  Insert a new entry into the table overwriting any existing match.
*/
extern int hdbm_replace(hdbm, datum, datum);
    
/*
  return true/false if key is present

  false is also returned if there are other errors
*/
extern int hdbm_probe(hdbm, datum);

/*
  return size of value if key is present
  
  if key is not present or there are other errors return -1
*/
extern int hdbm_size(hdbm, datum, int *);

/*
  Read the value associated with the key

  On success return datum in value and true

  On value not present or another error return false.
   
  The pointer returned in value points to a copy of the data which
  should be released with datum_free().
*/
extern int hdbm_read(hdbm, datum, datum *);

/*
  Delete the entry associated with key

  return

     1 if key was matched and delete successful
     0 otherwise
*/
extern int hdbm_delete(hdbm, datum);

/*
  Extract the value associated with the key
  (differs from read in that entry is deleted).

  On success return the value and true

  Otherwise return false.
*/
extern int hdbm_extract(hdbm, datum, datum *);
    
/*
  Second routine of iterator through keys in the database.

  On success return next key and true

  Otherwise return false
*/
extern int hdbm_next_key(hdbm, datum *);

/*
  First routine of iterator through keys in the database.

  On success return first key and true

  Otherwise return false
*/
extern int hdbm_first_key(hdbm, datum *);

/*
  Print the contents of the hash table out with user defined functions
  for printing keys and values.  If a function pointer is null a
  default function is provided.

  Return 1 on sucess, 0 on failure;
*/
extern int hdbm_print_table(hdbm db, void (*)(datum), void (*)(datum));

/*
  Print info about the hashtable

  Return 1 on sucess, 0 on failure.
*/
extern int hdbm_print_stats(hdbm);

/*
  Print info about resource usage (#calls, memory, I/O)
*/
extern void hdbm_print_usage();

#endif
