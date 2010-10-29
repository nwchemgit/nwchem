/*
 $Id$
*/
#ifndef _GET_WORD_H_
#define _GET_WORD_H_
/* get_word.h -
   Author - Eric Bylaska

*/
#include	<stdio.h>
#define	NIL	((char *) EOF)

extern char	*get_word(FILE *stream);
extern int 	get_line();
extern int 	to_eoln();
extern int 	get_int();
extern int 	get_string();
extern int 	get_float();
extern int 	get_end();
extern int 	remove_blanks();


#endif
