/*
 $Id: debug.c,v 1.1 2001-11-29 16:51:21 bylaska Exp $
   get_word.c -
   Author - Eric Bylaska

*/
#include	<stdio.h>
#include	"debug.h"


static	int	debug;

int debug_print()
{
      return debug;
}

void set_debug_print(int i)
{
   debug = i;
}
