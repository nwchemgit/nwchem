/*
 $Id: util_system.c,v 1.3 1999-02-16 07:42:04 d3e129 Exp $
 */

#include <stdio.h>
#include <unistd.h>

extern int system(const char *);

#ifdef CRAY
#include <fortran.h>
#endif

typedef long Integer;		/*  FORTRAN integer */

#ifdef CRAY
int fortchar_to_string(_fcd, int, char *, const int);
#else
int fortchar_to_string(const char *, int, char *, const int);
#endif

void ga_error(const char *, long);

#ifdef CRAY
Integer UTIL_SYSTEM(_fcd input)
{
    int lin  = _fcdlen(input);
#else
Integer util_system_(const char *input, int lin)
{
#endif
    char in[1024];
    if (!fortchar_to_string(input, lin, in, sizeof(in)))
	ga_error("util_system: fortchar_to_string failed for in",0);

#ifdef CRAY
    return 1;			/* Does not work on the Cray? */
#else
    return system(in);
#endif
}

