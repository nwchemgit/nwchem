/*
 $Id: util_system.c,v 1.5 1999-11-16 20:51:03 edo Exp $
 */

#include <stdio.h>
#ifndef WIN32
#include <unistd.h>
#endif

extern int system(const char *);

#ifdef CRAY
#include <fortran.h>
#define FATR
#endif
#ifdef WIN32
#include "typesf2c.h"
#endif

typedef long Integer;		/*  FORTRAN integer */

#if defined(CRAY) || defined(USE_FCD)
int fortchar_to_string(_fcd, int, char *, const int);
#else
int fortchar_to_string(const char *, int, char *, const int);
#endif

void ga_error(const char *, long);

#if defined(CRAY) || defined(USE_FCD)
Integer FATR UTIL_SYSTEM(_fcd input)
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

