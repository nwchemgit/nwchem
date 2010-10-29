/*$Id$*/
extern int context_set(const char *);
extern char *context_get(void);
extern int context_rtdb_store(int);
extern int context_rtdb_load(int);
extern int context_push(const char *);
extern int context_pop(const char *);
extern int context_rtdb_match(int, const char *, int, char *);
extern int context_prefix(const char *, char *, int);

 
#if defined(CRAY) || defined(WIN32)
#include "rtdb.cray.h"
#endif

