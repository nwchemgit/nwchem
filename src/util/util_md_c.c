/*
$Id: util_md_c.c,v 1.1 2000-07-24 17:46:58 d3j191 Exp $
*/
#include <time.h>
#include <sys/types.h>
#include <string.h>
#if defined(NEED_LOC)
int *loc_(int var)
/* Return address of var */
{return var;}
#endif
#if defined(LINUX)
#include <stdlib.h>
float rand_(int seed)
/* Return random number in [0.0:1.0] */
{float rnd = ((float) rand() / (float) RAND_MAX); return rnd;}
#endif
