/*
 $Id: nw_inp_from_string.c,v 1.7 2000-07-27 16:27:59 bjohnson Exp $
*/
#include "global.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifndef WIN32
#include <unistd.h>
#endif

#if defined(CRAY_T3E) || defined(CRAY_T3D) || defined(CRAY)
#define nw_inp_from_file_ NW_INP_FROM_FILE
#endif
extern int nw_inp_from_file_(int *, char *, int);
extern int string_to_fortchar(char *, int, char *);

int nw_inp_from_string(int rtdb, const char *input)
{
    char *filename = "temp.nw";
    FILE *file;
    char fstring[255];
    int status;

    if (ga_nodeid_() == 0) {
      if (!(file = fopen(filename,"w"))) {
	ga_error("nw_inp_from_string: failed to open temp.nw\n",0);
      }

      if (fwrite(input, 1, strlen(input), file) != strlen(input)) {
	ga_error("nw_inp_from_string: failed to write to temp.nw\n",0);
      }
      if (fwrite("\n", 1, 1, file) != 1) {
	ga_error("nw_inp_from_string: failed to write to temp.nw\n",0);
      }

      (void) fclose(file);
    }

    if (!string_to_fortchar(fstring, sizeof(fstring), filename)) {
      ga_error("nw_inp_from_string: fstring is too small",0);
    }

    status = nw_inp_from_file_(&rtdb, fstring, sizeof fstring);

    if (ga_nodeid_() == 0) (void) unlink(filename);

    return status;
}
