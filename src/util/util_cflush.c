/*
 $Id$
 */
#include <stdio.h>
/* this flushes stderr and stdout
   temporary hack for system where fortran flsuh is broken*/
int util_cflush_()
{
  fflush(stdout);
  fflush(stderr);
}
