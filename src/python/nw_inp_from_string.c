#include "global.h"
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#if defined(CRAY_T3E) || defined(CRAY_T3D)
#define NWC_NW_INP_FROM_FILE NW_INP_FROM_FILE
#else
#define NWC_NW_INP_FROM_FILE nw_inp_from_file_
#endif
extern int NWC_NW_INP_FROM_FILE(int *, char *, int);
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

      (void) fclose(file);
    }

    if (!string_to_fortchar(fstring, sizeof(fstring), filename)) {
      ga_error("nw_inp_from_string: fstring is too small",0);
    }

    status = NWC_NW_INP_FROM_FILE(&rtdb, fstring, sizeof fstring);

    if (ga_nodeid_() == 0) (void) unlink(filename);

    return status;
}
