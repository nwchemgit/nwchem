#include "typesf2c.h"

/* routine to convert a fortran string to a C string: */
/* Fortran callable version of f2cstring in global directory */

void c_cnvt_(Integer * clen, Integer * flen, char * cstr, char * fstr)
{
    int flenl, clenl;

    flenl = *flen; clenl = *clen;
    /* remove trailing blanks from fstr */
    while (flenl-- && fstr[flenl] == ' ') ;

    /* the postdecrement above went one too far */
    flenl++;

    /* truncate fstr to cstr size */
    if (flenl >= clenl)
        flenl = clenl - 1;

    /* ensure that cstr is NUL-terminated */
    cstr[flenl] = '\0';

    /* copy fstr to cstr */
    while (flenl--)
        cstr[flenl] = fstr[flenl];

}
