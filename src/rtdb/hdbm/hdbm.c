/*
 $Id$
 */

#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <string.h>
#ifndef WIN32
#include <unistd.h>
#endif
#include "hdbm.h"


/*
 *	On Cray replace the STDIO based IO with Cray's
 *	Flexable File IO (FFIO).
 */
#ifdef USE_FFIO
#include "stdiof2ffio.h"
#endif


static unsigned int hash(const void *, int);

#define MIN(a,b) ((a) <= (b) ? (a) : (b))

#define MAXHDBM 8		/* Max. no. of indep. data-bases */
#define NBINS 1021		/* No. of bins in the in-core hashtable */

#define HASH(a) (hash(a.dptr, a.dsize) % NBINS)
				/* Hash function for in-core hashtable */

typedef struct entry {
    struct entry *next;		/* Points down, null at bottom */
    datum key;			/* Always fully defined */
    datum value;		/* Defined only if NOT disk resident */
    long val_ptr;		/* Points to value on disk if disk resident */
    long rec_ptr;		/* Points to header record on disk if ... */
} entry;

typedef struct fflentry {	/* File free list entry */
    int size;
    long ptr;
} fflentry;

typedef struct {		/* Record that is stored on disk */
    int key_size;
    int val_size;
    int active;
} file_entry;

static struct {
    entry **table;		/* Will point to array of bins */
    int *nentry;		/* Will point to count of bin entries */
    int active;			/* Boolean ... true if open */
    char *name;			/* Name of db */
    int cur_bin;			/* Current bin for iterator */
    entry *cur_entry;		/* Pointer into current bin for iterator */
    FILE *file;			/* Stream for disk resident data bases */
    fflentry *ffl;		/* Will point to file free list */
    int nffl;			/* No. of entries in file free list */
} hash_tables[MAXHDBM];

static char cookie[] = "hdbm v1.0"; /* First characters in file */

static struct {                 /* Used to track memory leaks */
    unsigned long nmax;
    unsigned long n;
    int allocs;
    int frees;
} hdbm_malloc_stats;

static struct {			/* Curiosity */
    int reads;
    int inserts;
    int replaces;
    int extracts;
    int deletes;
    int reuses;
} call_stats;

static struct {			/* Used to track I/O operations */
    int reads;
    int writes;
    int seeks;
    unsigned long nread;
    unsigned long nwrite;
} io_stats;

/* Maximum size of the circular list of free file space.
   Note that the file is assumed dense so that space is only reused
   if it is an exact fit.  This occurs often within NWChem. */

#define MAX_FILE_FREE_LIST 512

/* Next few routines are wrappers for standard operations 
   so that have a hook to gather statistics, cache, ... */

static void *hdbm_malloc(size_t n)
{
    hdbm_malloc_stats.allocs++;
    hdbm_malloc_stats.n += n;
    if (hdbm_malloc_stats.nmax < hdbm_malloc_stats.n)
	hdbm_malloc_stats.nmax = hdbm_malloc_stats.n;
    return malloc(n);
}

static void hdbm_free(void *p, size_t n)
{
    hdbm_malloc_stats.frees++;
    hdbm_malloc_stats.n -= n;
    free(p);
}

static size_t hdbm_fread(void *pointer,
			 size_t size,
			 size_t num_items,
			 FILE *stream)
{
    io_stats.reads++;
    io_stats.nread += size*num_items;

    return fread(pointer, size, num_items, stream);
}

static size_t hdbm_fwrite(const void *pointer,
			  size_t size,
			  size_t num_items,
			  FILE *stream)
{
    io_stats.writes++;
    io_stats.nwrite += size*num_items;
    
    return fwrite(pointer, size, num_items, stream);
}


static int hdbm_fseek(FILE *stream,
		      long int offset,
		      int whence)
{
    int result;

    io_stats.seeks++;
    
    result = fseek(stream, offset, whence);
if ( result < 0 ) perror("hdbm_fseek");
    return ( result );
}

int hdbm_check(hdbm db)
/*
  Determine if db is a valid handle for an open hashtable.
  */
{
    return (db >= 0 && 
	    db < MAXHDBM && 
	    hash_tables[db].active);
}

void datum_wrap(const void *dptr, int dsize, datum *d)
/*
  Make a datum that refers to the region pointed to
  by the arguments .. NO MEMORY IS ALLOCATED.
  */
{
    d->dptr = (void *) dptr; d->dsize = dsize;
}

int datum_make(const void *dptr, int dsize, datum *d) 
/*
  Make a datum by duplicating the data described
  by the arguments. 
  
  Return  1 if successful
  0 if failed to allocate necessary memory
  */
{
    d->dsize = 0; d->dptr = 0;
    
    if (dsize && dptr) {
	if ((d->dptr = hdbm_malloc((size_t) dsize))) {
	    memcpy(d->dptr, dptr, (size_t) dsize);
	    d->dsize = dsize;
	    return 1;
	}
	else
	    return 0;
    }
    else
	return 1;			/* Keep quiet about empty datum */
}

int datum_duplicate(datum d, datum *r)
/*
  Duplicate the datum d, allocating space as required.  Return 1
  on success, 0 on failure.
  */
{
    void *dptr = d.dptr;
    int dsize  = d.dsize;
    
    r->dptr = 0; r->dsize = dsize; /* Note that datum can have non-zero
				      size but null pointer */
    
    if (dsize && dptr) {
	if ((r->dptr = hdbm_malloc((size_t) dsize))) {
	    memcpy(r->dptr, dptr, (size_t) dsize);
	    r->dsize = dsize;
	    return 1;
	}
	else
	    return 0;
    }
    else
	return 1;			/* Keep quiet about empty datum */
}

void datum_print_hex(datum d)
/*
  Print summary of datum to stdout in hex
  */
{
    int len = MIN(8, d.dsize);
    int i;
    unsigned char *ptr = d.dptr;
    
    (void) printf("(%d, %p, 0x",d.dsize, d.dptr);
    if (d.dptr) {
	for (i=0; i<len; i++) 
	    printf("%2.2x", ptr[i]);
	if (len < d.dsize)
	    printf(" ...");
    }
    printf(")");
    (void) fflush(stdout);
}

void datum_print_string(datum d)
/*
  Prints datum to stdout assuming a null terminated string.
  */
{
    int len = MIN(8, d.dsize-1);
    
    (void) printf("(%d, %p, \"",d.dsize, d.dptr);
    if (d.dptr) {
	(void) printf("%.*s", len, (char *) d.dptr);
	if (len < (d.dsize-1))
	    (void) printf("...");
    }
    (void) printf("\")");
    (void) fflush(stdout);
}

void datum_free(datum d)
/*
  Free the memory pointed to by d.dptr ... must have come from malloc.
  */
{
    if (d.dptr)
	hdbm_free(d.dptr, d.dsize);
}

static entry *NewEntry(datum key, datum value)
/*
  Assemble key, value into a new entry.  The datums (key and value)
  are duplicated
  */
{
    entry *new_entry;
    
    if (!(new_entry = (entry *) hdbm_malloc(sizeof(entry)))) {
	(void) fprintf(stderr, "HdbmNewEntry: malloc of entry failed\n");
	return (entry *) 0;
    }
    
    new_entry->next = (entry *) 0;
    if (!datum_duplicate(key, &new_entry->key))
	return (entry *) 0;
    if (!datum_duplicate(value, &new_entry->value))
	return (entry *) 0;
    
    return new_entry;
}

static void put_in_bin(hdbm db, entry *e) 
{
    int index = HASH(e->key);
    e->next = hash_tables[db].table[index];
    hash_tables[db].table[index] = e;
    hash_tables[db].nentry[index]++;
}

static int datum_fread(FILE *file, size_t file_ptr, size_t size, datum *d)
/*
  Allocate memory and read data from file returning result in datum.
  Return 1 on success, 0 on failure.
  */
{
    d->dsize = 0;
    
    if (!(d->dptr = hdbm_malloc((size_t) size))) {
	fprintf(stderr, "datum_fread: failed to malloc %lu\n", size);
	return 0;
    }
    if (hdbm_fseek(file, file_ptr, SEEK_SET)) {
	fprintf(stderr, "datum_fread: failed to position file %lu\n", 
		(long) file_ptr);
	return 0;
    }
    if (hdbm_fread((char *) d->dptr, (size_t) 1, (size_t) size, file) != size) {
	fprintf(stderr, "datum_fread: failed to read data %lu\n", size);
	return 0;
    }
    
    d->dsize = size;
    return 1;
}

static int hdbm_load(hdbm db, FILE *file)
/*
  Check that the file is indeed a data base and load data in
  */
{
    char header[32];
    file_entry fe;
    long rec_ptr;
    rewind(file);
    if ((hdbm_fread(header, sizeof(cookie), (size_t) 1, file) != 1) ||
	(strncmp(header, cookie, strlen(cookie)) != 0)) {
	(void) fprintf(stderr, "hdbm_load: cookie missing ... not a database\n");
	return 0;
    }

    /* Now simply keep reading until we hit end of file */
    
    rec_ptr = ftell(file);
    while (!hdbm_fseek(file, rec_ptr, SEEK_SET) &&
	   (hdbm_fread((char *) &fe, sizeof(fe), (size_t) 1, file) == 1)) {

	datum key, value;
	entry *e;

	if (fe.active) {
	    if (!datum_fread(file, rec_ptr+sizeof(fe), fe.key_size, &key)) {
		(void) fprintf(stderr, 
			       "hdbm_load: header present but key missing?\n");
		return 0;
	    }
	    
	    value.dptr = 0; value.dsize = fe.val_size;
	    
	    if (!(e = NewEntry(key, value)))
		return 0;
	    
	    datum_free(key); datum_free(value);
	    
	    e->rec_ptr = rec_ptr;
	    e->val_ptr = rec_ptr + ((long) sizeof(fe)) + fe.key_size;
	    
	    put_in_bin(db, e);
	}
	
	rec_ptr = rec_ptr + ((long) sizeof(fe)) + fe.key_size + fe.val_size;
	
    }
    
    /* Now check that we really are at EOF */
    
    (void) hdbm_fseek(file, 0L, SEEK_END);
    
    if (rec_ptr != ftell(file)) {
	fprintf(stderr, "hdbm_load: inconsistent end to file %ld != %ld\n",
		rec_ptr, ftell(file));
    }
    
    return TRUE;
}

int hdbm_open(const char *name, int use_file, hdbm *db)
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
{
    int i;
    int retval;
    
    *db = -1;			/* Invalid value */

    for (i=0; i<MAXHDBM; i++)
	if (hash_tables[i].active)
	    if (strcmp(name, hash_tables[i].name) == 0) {
		(void) fprintf(stderr,"hdbm_open: %s is already open\n",name);
		return 0;
	    }
    
    for (i=0; i<MAXHDBM; i++)	/* Find next free slot */
	if (!hash_tables[i].active)
	    break;
    if (i >= MAXHDBM) {
	(void) fprintf(stderr, 
		       "hdbm_open: too many hdbms have been opened (max=%d)\n",
		       MAXHDBM);
	return 0;
    }

    hash_tables[i].name = strdup(name);
    hash_tables[i].active = FALSE; /* So that errors return inactive db */
    
    if (!(hash_tables[i].table = 
	  (entry **) hdbm_malloc((size_t) (NBINS*sizeof(entry *))))) {
	(void) fprintf(stderr,"hdbm_open: malloc of table failed\n");
	return 0;
    }
    
    if (!(hash_tables[i].nentry = 
	  (int *) hdbm_malloc((size_t) (NBINS*sizeof(int))))){
	(void) fprintf(stderr, "hdbm_open: malloc of nentry failed\n");
	return 0;
    }
    
    (void) memset((void *) hash_tables[i].table, (char) 0, 
		  (size_t) (NBINS*sizeof(entry *)));
    (void) memset((void *) hash_tables[i].nentry, (char) 0,
		  (size_t) (NBINS*sizeof(int)));
    
    /* Now attempt to open the file if requested.
       "r+" will open an existing file for update.
       "w+" will create a new file or truncate an existing file. */
    
    if (use_file) {
	if (!(hash_tables[i].file = fopen(name, "r+b")))
	    if (!(hash_tables[i].file = fopen(name, "w+b"))) {
perror(name);
		(void) fprintf(stderr, "hdbm_open: open of %s failed\n", name);
		return 0;
	    }
	retval = hdbm_fseek(hash_tables[i].file, 0L, SEEK_END);
        if ( retval ) {
	    (void) fprintf(stderr, "hdbm_open: fseek on %s failed(%d)\n", name,retval);
	    return 0;
	}
	
	if (!(hash_tables[i].ffl = 
	      (fflentry *)hdbm_malloc(
		  (size_t)(MAX_FILE_FREE_LIST*sizeof(fflentry))))){
	    (void) fprintf(stderr, "hdbm_open: malloc of free list failed\n");
	    return 0;
	}
    
	hash_tables[i].nffl = 0;

	if (ftell(hash_tables[i].file) > 0) {
	    /* There is information in the file ... attempt to load it */
	    
	    if (!hdbm_load(i, hash_tables[i].file)) {
		hash_tables[i].active = FALSE;
		(void) fprintf(stderr,
			       "hdbm_open: failed to load file %s\n", name);
		return 0;
	    }
	}
	else {
	    /* The file is empty and the database is presumably new ...
	       put the cookie at the beginning */
	    
	    if (hdbm_fwrite(cookie, sizeof(cookie), 
			    (size_t) 1, hash_tables[i].file) != 1) {
		fprintf(stderr,
			"hdbm_open: failed to write cookie to %s\n", name);
		return 0;
	    }
	}
    }
    else {
	/* No file is to be used */
	hash_tables[i].file = (FILE *) 0;
    }
    
    hash_tables[i].active = TRUE;
    *db = i;
    return 1;
}

static int hdbm_clear(hdbm db)
/*
  Trash the contents of the db restoring to empty state.
  
  Return 1 on success, 0 on failure.
  */
{
    int bin, available;
    datum key;
    
    if (!hdbm_check(db)) {
	(void) fprintf(stderr,"hdbm_clear: invalid db handle passed %d\n", db);
	return 0;
    }
    
    for (available = hdbm_first_key(db, &key);
	 available;
	 available = hdbm_next_key(db, &key))  {
	if (!hdbm_delete(db, key)) {
	    (void) fprintf(stderr, "hdbm_clear: failed to delete key?\n");
	    datum_free(key);
	    return 0;
	}
	datum_free(key);
    }
    
    /* Sanity check ... are the bins really empty ? */
    
    for (bin=0; bin<NBINS; bin++)
	if (hash_tables[db].nentry[bin]) {
	    (void) fprintf(stderr,"hdbm_clear: ghost entries ?\n");
	    return 0;
	}
    
    return 1;
}

int hdbm_file_compress(const char *filename) 
/*
  Eliminate dead space in the database on disk by copying to and from
  an intermediate file ... the database in filename must be closed.
  */
{
    int l = (int) strlen(filename) + 5;
    char *tmp = hdbm_malloc((size_t) l);
    
    if (!tmp) {
	(void) fprintf(stderr,"hdbm_file_compress: no memory\n");
	return 0;
    }
    
    strcpy(tmp,filename);
    tmp[l-5] = '.';  tmp[l-4] = 't';  tmp[l-3] = 'm';  tmp[l-2] = 'p';
    tmp[l-1] = 0;

    if (!hdbm_file_copy(filename, tmp)) {
	(void) fprintf(stderr,
		       "hdbm_file_compress: copy to %s failed\n", tmp);
	return 0;
    }
    
    if (!hdbm_file_copy(tmp, filename)) {
	(void) fprintf(stderr,
		       "hdbm_file_compress: copy from %s failed\n", tmp);
	(void) fprintf(stderr,
		       "hdbm_file_compress: copy manually to recover\n");
	return 0;
    }
    
    (void) unlink(tmp);

    hdbm_free(tmp, l);
    
    return 1;
}

int hdbm_close(hdbm db)
/*
  Close this db ... the contents of any assoicated file are preserved.
  */
{
    if (!hdbm_check(db)) {
	(void) fprintf(stderr, 
		       "hdbm_close: invalid db handle passed %d\n", db);
	return 0;
    }
    if (hash_tables[db].file) {
	hdbm_free(hash_tables[db].ffl, MAX_FILE_FREE_LIST*sizeof(fflentry));

	(void) fclose(hash_tables[db].file);
	hash_tables[db].file = (FILE *) 0; /* So that delete called from clear
					      only destroys incore data*/
    }
    hdbm_clear(db);
    hdbm_free(hash_tables[db].table, NBINS*sizeof(entry *));
    hdbm_free(hash_tables[db].nentry, NBINS*sizeof(int));
    free(hash_tables[db].name);	/* This came from strdup */
    
    hash_tables[db].active = 0;
    
    return 1;
}

static int delete_file_entry(hdbm db, entry *e)
/*
  Delete the corresponding entry on disk.  Record it in the file_free_list.
  */
{
    FILE *file = hash_tables[db].file;
    file_entry fe;
    long act_ptr = e->rec_ptr + 
	(long) (((char *) &fe.active) - ((char *) &fe));
    int false = 0;
    int ind;

    if (hdbm_fseek(file, act_ptr, SEEK_SET)) {
	fprintf(stderr, "delete_file_entry: failed to position file %ld\n",
		act_ptr);
	return 0;
    }
    
    if (hdbm_fwrite((char *) &false, sizeof(false), (size_t) 1, file) != 1) {
	fprintf(stderr, "delete_file_entry: failed to inactivate header\n");
	return 0;
    }

    /* If room, put the file pointer into the free list for reuse */

    ind = hash_tables[db].nffl; /* No. of free list entries */
    if (ind < MAX_FILE_FREE_LIST) {
	hash_tables[db].ffl[ind].size = e->key.dsize+e->value.dsize;
	hash_tables[db].ffl[ind].ptr  = e->rec_ptr;                 
	hash_tables[db].nffl++;
	/* printf("saving %ld\n", e->rec_ptr); */
    }
    else {
	int i;	/* List is full ... if there is a smaller entry replace it */

	for (i=0; i<ind; i++) {
	    fflentry *t = hash_tables[db].ffl+i;
	    if ((e->value.dsize+e->key.dsize) > t->size) {
		t->size = e->key.dsize+e->value.dsize;
		t->ptr  = e->rec_ptr;                 
		/* printf("saving %ld\n", e->rec_ptr); */
		break;
	    }
	}
    }

    return 1;
}

static int write_file_entry(hdbm db, entry *e)
/*
  Store the key/value pair in the entry in free space of exactly the same 
  size or at the end of the file
  */
{
    FILE *file = hash_tables[db].file;
    file_entry fe;
    int ind, found=0;
    
    fe.key_size = e->key.dsize;
    fe.val_size = e->value.dsize;
    fe.active   = 1;

    /* Look for free space of the same size */

    for (ind=0; ind<hash_tables[db].nffl; ind++) {
	fflentry *tmp = hash_tables[db].ffl+ind;
	if (tmp->size == (fe.key_size +  fe.val_size)) {
	    if (hdbm_fseek(file, tmp->ptr, SEEK_SET)) {
		fprintf(stderr, "write_file_entry: seek failed%ld\n",
			tmp->ptr);
		return 0;
	    }

	    call_stats.reuses++;
	    found = 1;
	    hash_tables[db].nffl--;

	    /* Move last element over one just used to keep list dense */

	    *tmp = hash_tables[db].ffl[hash_tables[db].nffl];

	    break;
	}
    }

    if (!found) { /* No free space ... must append */
	/* printf("appending\n"); */
	if (hdbm_fseek(file, 0L, SEEK_END)) {
	    fprintf(stderr, "write_file_entry: failed to position file\n");
	    return 0;
	}
    }

    e->rec_ptr = ftell(file);
    
    if (hdbm_fwrite((char *) &fe, sizeof(fe), (size_t) 1, file) != 1) {
	fprintf(stderr, "write_file_entry: failed to write header\n");
	return 0;
    }
    
    if (hdbm_fwrite((char *) e->key.dptr, (size_t) 1, 
		    (size_t) e->key.dsize, file) != e->key.dsize) {
	fprintf(stderr, "write_file_entry: failed to write key\n");
	return 0;
    }
    
    e->val_ptr = ftell(file);
    
    if (hdbm_fwrite((char *) e->value.dptr, (size_t) 1, 
		    (size_t) e->value.dsize, file) != e->value.dsize) {
	fprintf(stderr, "write_file_entry: failed to write value\n");
	return 0;
    }
    
    return 1;
}

int hdbm_insert(hdbm db, datum key, datum value)
/*
  Insert a new entry into the table without looking
  for duplicates.  The datums (key and value) are duplicated.
  
  Return 1 on sucess, 0 on failure.
  
  */
{
    entry *top = NewEntry(key, value);
    
    call_stats.inserts++;
    
    if (!top)
	return 0;
    
    if (!hdbm_check(db)) {
	fprintf(stderr, "hdbm_insert: invalid database handle %d\n", db);
	return 0;
    }
    
    if (hash_tables[db].file) {
	/* If using a file write key/pair to file and free memory for value */
	
	int status = write_file_entry(db, top);
	datum_free(top->value); top->value.dptr = 0;
	if (!status) return 0;
    }
    
    put_in_bin(db, top);
    
    return 1;
}

int hdbm_replace(hdbm db, datum  key, datum value)
/*
  Insert a new entry into the table overwriting any existing match.
  
  Implemented as delete+insert so that only hdbm_insert ever
  inserts into the table.
  */
{
    call_stats.replaces++;
    
    if (!hdbm_check(db)) {
	fprintf(stderr, "hdbm_replace: invalid database handle %d\n", db);
	return 0;
    }
    
    (void) hdbm_delete(db, key);
    
    return hdbm_insert(db, key, value);
}

static void search(entry **table, datum key, int *p_index, 
		   entry **p_cur, entry **p_prev)
/*
  Search for a match for key in the referenced table
  
  On success
  
  index  returns the bucket no.
  p_cur  points to the desired entry
  p_prev points to previous entry in linked list
  of the bucket
  
  On failure
  
  p_cur and p_prev return 0
  index returns junk
  */
{
    int index = HASH(key);
    entry *cur, *prev;
    
    *p_index = index;
    
    for (prev=0, cur=table[index]; cur; prev=cur, cur=cur->next)
	if (DATUMS_MATCH(key, cur->key)) {
	    *p_cur = cur;
	    *p_prev = prev;
	    return;
	}
    
    *p_cur = *p_prev = (entry *) 0;
}

int hdbm_probe(hdbm db, datum key)
/*
  return true/false if key is present
  
  false is also returned if there are other errors
  */
{
    int index;
    entry *cur, *prev;
    
    if (!hdbm_check(db)) {
	fprintf(stderr, "hdbm_probe: invalid database handle %d\n", db);
	return 0;
    }
    
    search(hash_tables[db].table, key, &index, &cur, &prev);
    
    return (int) (cur != 0);
}

int hdbm_size(hdbm db, datum key, int *size)
/*
  Return 1 and size of value if key is present.
  
  Return 0 otherwise.
  */
{
    int index;
    entry *cur, *prev;
    
    *size = 0;
    
    if (!hdbm_check(db)) {
	fprintf(stderr, "hdbm_size: invalid database handle %d\n", db);
	return 0;
    }
    
    search(hash_tables[db].table, key, &index, &cur, &prev);
    
    if (cur) {
	*size =  cur->value.dsize;
	return 1;
    }
    else
	return 0;
}

int hdbm_read(hdbm db, datum key, datum *value)
/*
  Read the value associated with the key
  
  On success return datum in value and true
  
  On value not present or another error return false.
  
  The pointer returned in value points to a copy of the data which
  should be released with datum_free().
  */
{
    int index;
    entry *cur, *prev;
    
    call_stats.reads++;
    
    value->dsize = 0; value->dptr = 0;
    
    if (!hdbm_check(db)) {
	fprintf(stderr, "hdbm_read: invalid database handle %d\n", db);
	return 0;
    }
    
    search(hash_tables[db].table, key, &index, &cur, &prev);
    
    if (cur) {
	if (hash_tables[db].file)
	    return datum_fread(hash_tables[db].file, 
			       cur->val_ptr, cur->value.dsize, value);
	else
	    return datum_duplicate(cur->value, value);
    }
    else 
	return 0;
}


int hdbm_file_copy(const char *filefrom, const char *fileto)
/*
  Delete any original content of fileto and then copy the 
  database in filefrom into fileto.  Both files are closed
  at the end of the operation.
  */
{
    hdbm dbfrom, dbto;
    datum key, value;
    int available;
    FILE *file;

    /* Truncate the to-file to zero length or create it ... this is
       better than unlinking it since it may be a link or a special file */
    
    if (!(file = fopen(fileto, "w"))) {
	(void) fprintf(stderr,"hdbm_file_copy: failed to truncate (to) %s\n", 
		       fileto);
	return 0;
    }
    (void) fclose(file);
    
    /* Open the databases */

    if (!hdbm_open(filefrom, 1, &dbfrom)) {
	(void) fprintf(stderr,"hdbm_file_copy: failed to open (from) %s\n", 
		       filefrom);
	return 0;
    }
    if (!hdbm_open(fileto, 1, &dbto)) {
	(void) fprintf(stderr,"hdbm_file_copy: failed to open (to) %s\n", fileto);
	(void) hdbm_close(dbfrom);
	return 0;
    }
    
    /* Do the copy */
    
    for (available = hdbm_first_key(dbfrom, &key);
	 available;
	 available = hdbm_next_key(dbfrom, &key)) {
	
	int ok = hdbm_read(dbfrom, key, &value);
	if (!ok)
	    (void) fprintf(stderr, "hdbm_file_copy: failed reading from %s\n", 
			   filefrom);
	else {
	    ok = hdbm_insert(dbto, key, value);
	    datum_free(key); datum_free(value);
	    if (!ok)
		(void) fprintf(stderr, "hdbm_file_copy: failed writing to %s\n", 
			       fileto);
	}
	if (!ok) {
	    (void) hdbm_close(dbfrom);
	    (void) hdbm_close(dbto);
	    return 0;
	}
    } 
    
    return (hdbm_close(dbfrom) && hdbm_close(dbto));
}

int hdbm_delete(hdbm db, datum key)
/*
  Delete the entry associated with key
  
  return
  
  1 if key was matched and delete successful
  0 otherwise
  */
{
    int index;
    entry *cur, *prev;
    
    call_stats.deletes++;
    
    if (!hdbm_check(db)) {
	fprintf(stderr, "hdbm_delete: invalid database handle %d\n", db);
	return 0;
    }
    
    search(hash_tables[db].table, key, &index, &cur, &prev);
    
    if (cur) {
	if (prev)
	    prev->next = cur->next;
	else
	    hash_tables[db].table[index] = cur->next;
	datum_free(cur->value);
	
	if (hash_tables[db].file)
	    delete_file_entry(db, cur);
	
	datum_free(cur->key);	/* Free up memory for key/value/entry */
	datum_free(cur->value);
	hdbm_free((void *) cur, sizeof(*cur));
	
	hash_tables[db].nentry[index]--;
	return 1;
    }
    else
	return 0;
}

int hdbm_extract(hdbm db, datum key, datum *value)
/*
  Extract the value associated with the key
  (differs from read in that entry is deleted).
  
  On success return the value and true
  
  Otherwise return false.
  
  Implemented as read+delete so that only delete/insert manipulate entries.
  */
{
    call_stats.extracts++;
    
    if (!hdbm_check(db)) {
	fprintf(stderr, "hdbm_extract: invalid database handle %d\n", db);
	return 0;
    }
    
    value->dsize = 0; value->dptr = 0;
    
    if (!hdbm_read(db, key, value))
	return 0;
    (void) hdbm_delete(db, key);
    return 1;
}

static unsigned int hash(const void *vdata, int n)
{
/*  unsigned int multiplier = 257; */
    unsigned int multiplier = 1033;
    unsigned int value = 0;
    const unsigned char *data = vdata;
    
    while (n--)
	value = *data++ + multiplier * value;
    
    return value;
}

int hdbm_next_key(hdbm db, datum *key)
/*
  Second routine of iterator through keys in the database.
  
  On success return next key and true
  
  Otherwise return false
  */
{
    entry *cur_entry;
    int cur_bin;
    
    key->dptr = 0; key->dsize = 0;
    
    if (!hdbm_check(db)) {
	fprintf(stderr,"hdbm_next_key: invalid database handle %d\n", db);
	return 0;
    }
    
    cur_entry = hash_tables[db].cur_entry;
    cur_bin = hash_tables[db].cur_bin;
    
    if ((cur_bin < 0) || (cur_bin >= NBINS))
	return 0;
    
    if (!cur_entry) {
	while (cur_bin < (NBINS-1)) {
	    cur_bin++;
	    if ((cur_entry = hash_tables[db].table[cur_bin]))
		break;
	}
    }
    
    hash_tables[db].cur_bin = cur_bin;
    
    if (cur_entry) {
	hash_tables[db].cur_entry = cur_entry->next;
	return datum_duplicate(cur_entry->key, key);
    }
    else {
	hash_tables[db].cur_entry = 0;
	return 0;
    }
}

int hdbm_first_key(hdbm db, datum *key)
/*
  First routine of iterator through keys in the database.
  
  On success return first key and true
  
  Otherwise return false
  */
{
    entry *cur_entry = 0;
    int cur_bin;
    
    key->dptr = 0; key->dsize = 0;
    
    if (!hdbm_check(db)) {
	fprintf(stderr, "hdbm_first_key: invalid database handle %d\n", db);
	return 0;
    }
    
    for (cur_bin=0; 
	 (cur_bin < NBINS) && (hash_tables[db].nentry[cur_bin] == 0); 
	 cur_bin++)
	;
    
    hash_tables[db].cur_bin = cur_bin;
    
    if (cur_bin != NBINS) {
	cur_entry = hash_tables[db].table[cur_bin];
	hash_tables[db].cur_entry = cur_entry->next;
	return datum_duplicate(cur_entry->key, key);
    }
    else {
	hash_tables[db].cur_entry = (entry *) 0;
	return 0;
    }
}

int hdbm_print_table(hdbm db, 
		     void (*print_key)(datum), void (*print_value)(datum))
/*
  Print the contents of the hash table out with user defined functions
  for printing keys and values.  If a function pointer is null a
  default function is provided.
  
  Return 1 on sucess, 0 on failure;
  */
{
    datum key, value;
    int n, available;
    
    if (!hdbm_check(db)) {
	fprintf(stderr, "hprint_table: invalid database handle %d\n", db);
	return 0;
    }
    
    printf("Contents of \"%s\"\n\n",hash_tables[db].name);
    
    /* Bind default functionality ... */
    
    if (!print_key) print_key = datum_print_hex;
    if (!print_value) print_value = datum_print_hex;
    
    n = 0;
    available = hdbm_first_key(db, &key);
    
    while (available) {
	printf("%6d key=", n); print_key(key); 
	if (!hdbm_read(db, key, &value)) {
	    (void) fprintf(stderr, "hdbm_print_table: got key but no value?\n");
	    return 0;
	}
	printf(", value="); print_value(value); printf("\n");
	datum_free(key); datum_free(value);
	n++;
	available = hdbm_next_key(db, &key);
    }
    
    if (n == 0)
	printf("Empty table\n");
    
    printf("\n");
    
    return 1;
}    

void hdbm_print_usage(void)
{
    printf("\n HDBM usage statistics\n ---------------------\n");
    printf(" Reads         = %d\n", call_stats.reads);
    printf(" Inserts       = %d\n", call_stats.inserts);
    printf(" Replaces      = %d\n", call_stats.replaces);
    printf(" Extracts      = %d\n", call_stats.extracts);
    printf(" Deletes       = %d\n", call_stats.deletes);
    printf(" Reuses        = %d\n", call_stats.reuses);
    printf("\n");
    printf(" Memory allocs = %d\n", hdbm_malloc_stats.allocs);
    printf(" Memory frees  = %d\n", hdbm_malloc_stats.frees);
    printf(" Memory max    = %lu\n", hdbm_malloc_stats.nmax);
    printf(" Memory in use = %lu\n", hdbm_malloc_stats.n);
    printf("\n");
    printf(" IO seeks      = %d\n", io_stats.seeks);
    printf(" IO reads      = %d\n", io_stats.reads);
    printf(" IO writes     = %d\n", io_stats.writes);
    printf(" IO nread      = %lu\n",io_stats.nread);
    printf(" IO nwrite     = %lu\n",io_stats.nwrite);
    printf("\n");
    
    fflush(stdout);
}

int hdbm_print_stats(hdbm db)
/*
  Print info about the hashtable
  
  Return 1 on sucess, 0 on failure.
  */
{
    int bin, nprint;
    
    if (!hdbm_check(db)) {
	fprintf(stderr, "hdbm_print_stats: invalid database handle %d\n", db);
	return 0;
    }
    
    printf(" Hash table bin statistics for \"%s\"\n\n",hash_tables[db].name);
    
    for (bin=0, nprint=0; bin<NBINS; bin++) {
	nprint += printf(" %3d", hash_tables[db].nentry[bin]);
	if (nprint > 72) {
	    printf("\n");
	    nprint = 0;
	}
    }
    if (nprint) printf("\n");
    
    return 1;
}    

int hdbm_file_flush(hdbm db)
/*
  Flush any pending I/O to disk
  */
{
    if (!hdbm_check(db)) {
	fprintf(stderr, "hdbm_flush_file: invalid database handle %d\n", db);
	return 0;
  }
    
    if (hash_tables[db].file) {
#ifdef USE_FFIO
	if (stdiof_fflush(hash_tables[db].file)) {
#else
	if (fflush(hash_tables[db].file)) {
#endif
	    (void) fprintf(stderr,"hdbm_file_flush: failed to flush %s\n",
			   hash_tables[db].name);
	    (void) fflush(stderr);
	    return 0;
	}
    }
    else {
	(void) fprintf(stderr,"hdbm_file_flush: no file for %s\n",
		       hash_tables[db].name);
	(void) fflush(stderr);
	return 0;
    }
    
    return 1;
}

