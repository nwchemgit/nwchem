/*
 $Id: push_inp_cstring.c 19695 2013-03-08 16:51:02Z d3y133 $
*/
#include "ga.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef CRAY_T3E
#define FATR
#include <fortran.h> /* Required for Fortran-C string interface on Crays */
#endif
#ifndef WIN32
#include <unistd.h>
#else
#include "typesf2c.h"
#endif

#if defined(CRAY_T3E)  || defined(WIN32)
#define push_inp_string_ PUSH_INP_STRING
#endif

#if defined(CRAY_T3E) || defined(USE_FCD) || defined(WIN32)
extern Integer FATR push_inp_string_(_fcd inp);
#else
extern Integer FATR push_inp_string_(char *inp, int len);
#endif

int push_inp_cstring(const char *input)
{
    int status;
#if defined(USE_FCD) || defined(CRAY_T3E) || defined(WIN32)
    _fcd inp;
#else
    const char *inp = input;
#endif

#if defined(CRAY_T3E)
    inp = _cptofcd(input, strlen(input));
    status = push_inp_string_(inp);
#elif defined(WIN32)
    inp.string = input;
    inp.len = strlen(input);
    status = push_inp_string_(inp);
#elif defined(USE_FCD)
#error Do something about _fcd
#else
    status = push_inp_string_(inp, strlen(inp));
#endif
    return status;
}
