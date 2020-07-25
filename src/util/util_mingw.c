/*  $Id$*/
/* stubs for functions missing in MinGW */
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
void bzero(void *ptr, size_t n)
{
  memset(ptr, 0, n);
}
void nice(int n)
{
/* do nothing */
}
void bcopy(const void *src, void *dest, size_t n)
{
memcpy(dest, src, n);
}
int pause(void)
{
  fprintf(stderr, " pause() function not available \n");
  return -1;
}
