/*
 $Id: name.c,v 1.1 2001-08-30 16:58:36 bylaska Exp $
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
   else s = "?";
   return s;
}

