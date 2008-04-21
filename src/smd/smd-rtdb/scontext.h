/*$Id: scontext.h,v 1.1 2008-04-21 18:19:27 marat Exp $*/
extern int context_set(const char *);
extern char *context_get(void);
extern int context_srtdb_store(int);
extern int context_srtdb_load(int);
extern int context_push(const char *);
extern int context_pop(const char *);
extern int context_srtdb_match(int, const char *, int, char *);
extern int context_prefix(const char *, char *, int);

 
#if defined(CRAY) || defined(WIN32)
#include "srtdb.cray.h"
#endif

