/*
 $Id: c_inp.c,v 1.6 1999-11-13 02:22:25 bjohnson Exp $
 */

#ifdef CRAY
#include <fortran.h>
#endif
#include "typesf2c.h"
/* routine to convert a fortran string to a C string: */
/* Fortran callable version of f2cstring in global directory */

#if defined(CRAY) || defined(USE_FCD)
void FATR C_CNVT(clen, flen, cfcd, ffcd)
Integer *clen, *flen;
_fcd cfcd;
_fcd ffcd;
{
    char *cstr = _fcdtocp(cfcd);
    char *fstr = _fcdtocp(ffcd);
#else
void c_cnvt_(clen, flen, cstr, fstr)
Integer *clen, *flen;
char *cstr;
char *fstr;
{
#endif
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
