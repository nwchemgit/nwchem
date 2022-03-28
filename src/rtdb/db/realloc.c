/*$Id$*/
#include <sys/types.h>

#include <stdlib.h>

void *
__fix_realloc(p, n)
	void *p;
	size_t n;
{
#if defined(IBM)
void *calloc(size_t,size_t), *realloc(void *, size_t);
#endif
        return ((p == (void *)0) ? calloc(n, (size_t) 1) : realloc(p, n));
}
