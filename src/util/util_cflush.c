/*
 $Id$
 */
#include <stdio.h>
/* this flushes stderr and stdout
   temporary hack for system where fortran flush is broken

   Also useful for GNU compilers where the Fortran I/O and
   C I/O are buffered independently. Hence one needs to 
   flush both I/O systems.
*/
int util_cflush_()
{
  fflush(stdout);
  fflush(stderr);
  return 0;
}
