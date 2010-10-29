/*
   $Id$
*/

#include        <stdlib.h>
#include        <stdio.h>
#include        <string.h>
#include        <math.h>

#include        "paw_sdir.h"



/*internal core data*/
static char sdir_name[500];
static int  debug;

void  paw_set_debug(int debug_in)
{
    debug = debug_in;
}
int paw_debug()
{
    return debug;
}

void  paw_set_sdir(char *sdir,int m9)

{
    (void) strncpy(sdir_name, sdir, m9);
    sdir_name[m9] = '\0';
    strcat(sdir_name,"/");
    sdir_name[m9+1] = '\0';
}

char *paw_sdir()
{
    return sdir_name;
}

