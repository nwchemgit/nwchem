#include <sys/types.h>
#include <sys/time.h>

#if defined(CRAY)
void CPUTM(ai)
#else
void cputm_(ai)
#endif
int *ai;
{
struct timeval tp;
struct timezone tzp;
int i;

	gettimeofday(&tp,&tzp);

	/* B
	printf("seconds time=%ld\n",tp.tv_sec);
	printf("microseconds time=%ld\n",tp.tv_usec);
	E */
	i = tp.tv_sec & 0xffffff;
	i = i*100;
	i += (tp.tv_usec / 10000);
	*ai = i;
}
