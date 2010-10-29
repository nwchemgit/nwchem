/*
 $Id$
*/
#include <stdio.h>
void sgi_flush6_(void)
{
  int return_code;
  return_code = fflush(stdout);
  if (return_code == 0) return;
  (void) perror("Error detected in sgi_flush6:");
}
