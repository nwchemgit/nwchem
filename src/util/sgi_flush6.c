/*
 $Id: sgi_flush6.c,v 1.2 1999-07-27 21:00:19 d3e129 Exp $
*/
#include <stdio.h>
void sgi_flush6_(void)
{
  int return_code;
  return_code = fflush(stdout);
  if (return_code == 0) return;
  (void) perror("Error detected in sgi_flush6:");
}
