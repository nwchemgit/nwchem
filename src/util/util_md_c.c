/*
$Id$
*/
#include <time.h>
#include <sys/types.h>
#include <string.h>
#if defined(NEED_LOC)
void *loc_(void* var)
/* This routine is called by Fortran which passes by address (ie. pointer).
   we then pass that back by value so that Fortran can learn the "location"
   of that variable in memory.  The size of the Integer than Fortran expects
   as the return value need to match the pointer's size  */
{return var;}
#endif
#if defined(LINUX)
#include <stdlib.h>
float rand_(int seed)
/* Return random number in [0.0:1.0] */
{float rnd = ((float) rand() / (float) RAND_MAX); return rnd;}
#endif
