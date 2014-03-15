/*
 $Id$
   get_word.c -
   Author - Eric Bylaska

*/
#include	<stdio.h>
#include	<string.h>
#include	"get_word.h"


static	char	word[72];

char	*get_word(FILE *stream)
{
    if (fscanf(stream,"%s",word) != EOF)
        return word;
    else
        return NIL;
}


int get_line(stream,line,maxlen)
FILE *stream;
char *line;
int  maxlen;
{
    int   c,i;

    for (i=0; i<(maxlen-1) && (c=fgetc(stream))!=EOF && c!='\n'; ++i)
        line[i] = c;

    if (c=='\n')
    {
        line[i] = c;
        ++i;
    }
    line[i] = '\0';
    return i;
}


int to_eoln(stream)
FILE *stream;
{
    int   c,i;

    i = 0;
    while ((c=fgetc(stream))!=EOF && c!='\n')
        ++i;

    return i;
}


int     get_int(FILE *stream, int *ii)
{
    int c,value,n,i;

    for (i=0; i<72; ++i) word[i] = '\0';
    value = remove_blanks(stream);

    value = 0;
    n     = 0;
    while ((c=fgetc(stream))!=EOF && c!='\n' && c!=' ' && c!='\t')
    {
        word[n] = c;
        ++n;
    }
    if ((c=='\n') || (c==' ') || (c=='\t'))
    {
        ungetc(c,stream);
        --n;
    }

    value = sscanf(word,"%d",ii);
    if (!value) for (i=0; i<n; ++i) ungetc(word[n-i],stream);


    return value;
}




int get_float(FILE *stream, double *ff)
{
    int c,value,n,i;

    for (i=0; i<72; ++i) word[i] = '\0';
    value = remove_blanks(stream);

    value = 0;
    n     = 0;
    while ((c=fgetc(stream))!=EOF && c!='\n' && c!=' ' && c!='\t')
    {
        word[n] = c;
        ++n;
    }
    if ((c=='\n') || (c==' ') || (c=='\t'))
    {
        ungetc(c,stream);
        --n;
    }

    value = sscanf(word,"%lf",ff);
    if (!value) for (i=0; i<n; ++i) ungetc(word[n-i],stream);


    return value;
}





int     get_string(FILE *stream, char *string)
{
    int c,value,n,i;

    for (i=0; i<72; ++i) word[i] = ' ';
    value = remove_blanks(stream);

    n     = 0;
    while ((c=fgetc(stream))!=EOF && c!='\n' && c!=' ' && c!='\t')
    {
        word[n] = c;
        ++n;
    }

    if ((c=='\n'))
    {
        ungetc(c,stream);
        --n;
    }

    value = sscanf(word,"%s",string);
    if (!value) for (i=0; i<n; ++i) ungetc(word[n-i],stream);


    return value;
}


int remove_blanks(FILE *stream)
{
    int c,n=0,value;

    while ((c=fgetc(stream))!=EOF && c!='\n' && (c==' ' || c=='\t'))
        ++n;

    ungetc(c,stream);


    value = 1;
    if ((c=='\n'))
    {
        value = 0;
        --n;
    }

    return value;
}


int     get_end(FILE *stream)
{
    int c,value,n,i;

    for (i=0; i<72; ++i) word[i] = '\0';
    value = remove_blanks(stream);

    value = 0;
    n     = 0;
    while ((c=fgetc(stream))!=EOF && c!='\n' && c!=' ' && c!='\t')
    {
        word[n] = c;
        ++n;
    }
    if ((c=='\n') || (c==' ') || (c=='\t'))
    {
        ungetc(c,stream);
        --n;
    }
    for (i=0; i<n; ++i) ungetc(word[n-i],stream);

    value = !strcmp("<end>",word);


    return value;
}

