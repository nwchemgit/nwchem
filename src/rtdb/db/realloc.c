/*$Id: realloc.c,v 1.1 1995-03-31 01:55:36 d3g681 Exp $*/
#include <sys/types.h>

#include <stdlib.h>

void *
__fix_realloc(p, n)
	void *p;
	size_t n;
{
#if defined(SGI) || defined(IBM)
void *calloc(size_t,size_t), *realloc(void *, size_t);
#endif
        return ((p == (void *)0) ? calloc(n, (size_t) 1) : realloc(p, n));
}
