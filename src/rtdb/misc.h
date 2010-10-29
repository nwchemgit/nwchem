/*$Id$*/
#define error(format,data) \
  {(void) fflush(stdout); \
   (void) fprintf(stderr,format,data); \
   (void) fprintf(stderr, "error detected at line %d of file %s\n", \
		  __LINE__, __FILE__); \
   (void) fflush(stderr); \
   (void) abort();}


/* Macro returns pointer of type (type *) to n_items objects of
   type (type).  Error is called with the given error message
   if malloc fails.  Malloc, rather than calloc, is used as this
   is much faster.  Thus the memory is not initialized. */

#ifdef ALLOCATE_PRINT

#define ALLOCATE(pointer, n_items, type, err_msg) \
    {(void) printf("ALLOCATE:file %s:line %d: %s: %d", \
		   __FILE__, __LINE__, err_msg,(n_items)*sizeof(type)); \
     (void) fflush(stdout); \
     pointer = (type *) malloc((unsigned) ((n_items)*sizeof(type))); \
     (void) printf(": %p\n",pointer); \
     (void) fflush(stdout); \
     if (!pointer && n_items) \
       error(err_msg, (int) (n_items)); }

#define FREE(pointer) \
   {(void) printf("FREE:file %s:line %d: %x\n", __FILE__, __LINE__, pointer); \
    if (pointer) \
      free((char *) pointer); \
    else \
      (void) printf("Freeing NULL pointer...\n");}

#else

#define ALLOCATE(pointer, n_items, type, err_msg) \
  { pointer = (type *) malloc((unsigned) ((n_items)*sizeof(type))); \
    if (!pointer && n_items) \
      error(err_msg, (int) n_items); }

#define FREE(pointer) \
    { if (pointer) \
      free((char *) pointer); \
    else \
      (void) printf("Freeing NULL pointer...\n") ; }

#endif

