/*$Id: context.h,v 1.3 1995-02-02 23:21:59 d3g681 Exp $*/
extern int context_set(const char *);
extern char *context_get(void);
extern int context_rtdb_store(int);
extern int context_rtdb_load(int);
extern int context_push(const char *);
extern int context_pop(const char *);
extern int context_rtdb_match(int, const char *, int, char *);
extern int context_prefix(const char *, char *, int);

 
#ifdef CRAY
#include "rtdb.cray.h"
#endif

