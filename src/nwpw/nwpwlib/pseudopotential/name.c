/*
 $Id$
*/
#include	"name.h"

char *spd_Name(int l)
{
    char *s;
    if (l==0)      s = "s";
    else if (l==1) s = "p";
    else if (l==2) s = "d";
    else if (l==3) s = "f";
    else if (l==4) s = "g";
    else if (l==5) s = "h";
    else if (l==6) s = "i";
    else if (l==7) s = "j";
    else s = "?";
    return s;
}

