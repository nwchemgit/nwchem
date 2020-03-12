/*
 $Id$
 */

#include <sys/types.h>
#include <stdio.h>
#ifdef WIN32
#include "typesf2c.h"
#else
#include <sys/time.h>
#endif


#if (defined(CRAY) || defined(WIN32))&& !defined(__crayx1) && !defined (__MINGW32__)
#ifndef WIN32
#define FATR
#endif
void FATR CPUTM(ai)
#else
void cputm_(ai)
#endif
int *ai;
{
#ifndef WIN32
  /* !!! Comment out function for WIN32 just so we can keep going */
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
#else
	printf("selci/cputm.c not implemented yet for WIN32\n");
#endif
}
