#include <stdio.h>
#include <string.h>

void f_memzero_(void* address, int* length)
{
  memset(address,0.0,*length);
  return;
}

