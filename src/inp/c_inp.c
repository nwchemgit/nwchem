#include "types.f2c.h"
/* routine to convert a fortran string to a C string: */
/* Fortran callable version of f2cstring in global directory */

c_cnvt_(clen, flen, cstr, fstr)
Integer *clen, *flen;
char *cstr;
char *fstr;
{
    /* remove trailing blanks from fstr */
    while ((*flen)-- && fstr[*flen] == ' ') ;

    /* the postdecrement above went one too far */
    (*flen)++;

    /* truncate fstr to cstr size */
    if (*flen >= *clen)
        *flen = *clen - 1;

    /* ensure that cstr is NUL-terminated */
    cstr[*flen] = '\0';

    /* copy fstr to cstr */
    while ((*flen)--)
        cstr[*flen] = fstr[*flen];

}
