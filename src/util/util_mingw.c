/*  $Id: win32_cpu.c 19707 2010-10-29 17:59:36Z d3y133 $*/
/* stubs for functions missing in MinGW */
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
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
