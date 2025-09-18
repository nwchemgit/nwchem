#include <stdio.h>
#ifndef WIN32 
#include <unistd.h>
#endif

extern int system(const char *);

#if defined(WIN32) && !defined(__MINGW32__)
#include "typesf2c.h"
#endif

void GA_Error(char *str, int code);

typedef long Integer;		/*  FORTRAN integer */

int fortchar_to_string(const char *, int, char *, const int);

Integer util_system_(const char *input, int lin)
{
#if defined(LINUX) && !defined(__x86_64__)
    int i;
    void (*Siginit)();
#endif
    char in[1024];
    if (!fortchar_to_string(input, lin, in, sizeof(in)))
	GA_Error("util_system: fortchar_to_string failed for in",0);

#if defined(CATAMOUNT)
    GA_Error("system calls do not work on this machine", 0);
#else
    return system(in);
#endif
}
