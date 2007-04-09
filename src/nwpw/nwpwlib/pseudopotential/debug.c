/*
 $Id: debug.c,v 1.2 2007-04-09 22:55:51 d3p708 Exp $
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
