/*
 $Id: util_cflush.c,v 1.1 2004-10-06 05:16:45 edo Exp $
 */
#include <stdio.h>
/* this flushes stderr and stdout
   temporary hack for system where fortran flsuh is broken*/
int util_cflush_()
{
  fflush(stdout);
  fflush(stderr);
}
