/*
 $Id: get_word.c,v 1.1 2001-08-30 16:58:35 bylaska Exp $
   get_word.c -
   Author - Eric Bylaska

*/
#include	<stdio.h>
#include	"get_word.h"


static	char	word[50];

char	*get_word(FILE *stream)
{
   if (fscanf(stream,"%s",word) != EOF)
      return word;
   else
      return NIL;
}

