/*
 $Id: nw_inp_from_string.c,v 1.9 2000-08-26 19:58:09 d3g681 Exp $
*/
#include "global.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifndef WIN32
#include <unistd.h>
#endif
#include "typesf2c.h"

#if defined(CRAY_T3E) || defined(CRAY_T3D) || defined(CRAY) || defined(WIN32)
#define nw_inp_from_file_ NW_INP_FROM_FILE
#endif

#if defined(CRAY) || defined(USE_FCD) || defined(WIN32)
extern Integer FATR nw_inp_from_file_(Integer *rtdb, _fcd filename);
#else
extern Integer FATR nw_inp_from_file_(Integer *rtdb, char *filename, int flen);
#endif

int nw_inp_from_string(Integer rtdb, const char *input)
{
    char *filename = "temp.nw";
    FILE *file;
#if defined(USE_FCD) || defined(CRAY) || defined(WIN32)
    _fcd fstring;
#else
    char fstring[255];
#endif
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

#if defined(CRAY)
      _fstring = _cptofcd(filename, strlen(filename));
#elif defined(WIN32)
    fstring.string = filename;
    fstring.len = strlen(filename);
    status = nw_inp_from_file_(&rtdb, fstring);
#elif defined(USE_FCD)
#error Do something about _fcd
#else
    status = nw_inp_from_file_(&rtdb, filename, strlen(filename));
#endif


    if (ga_nodeid_() == 0) (void) unlink(filename);

    return status;
}
