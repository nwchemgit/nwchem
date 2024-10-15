/*
  Try to route all sequential file operations thru eaf/elio
  so that can handle extents correctly
*/


#include <stdio.h>
#ifdef WIN32
#include <io.h>
#define F_OK 2
#include "typesf2c.h"
#else
#include <unistd.h>
#endif

#include "eaf.h"
#include "ga.h"

int fortchar_to_string(const char *, int, char *, const int);

void util_file_unlink(const char *filename)
/*
  Delete the file.  If the file does not exist, quietly return.
  If the file exists and the unlink fails then abort.
  */
{
  /*
    if (access(filename, F_OK) == 0) {
	if (unlink(filename)) {
	    fprintf(stderr,"util_file_unlink: failed unlinking %s\n",
		    filename);
	    GA_Error("util_file_unlink",0);
	}
    }
  */
  if (EAF_Delete(filename) != 0)
    GA_Error("util_file_unlink",0);
}

void util_file_unlink_(const char *input, int lin)
{
    char in[255];
    if (!fortchar_to_string(input, lin, in, sizeof(in)))
	GA_Error("util_file_unlink: fortchar_to_string failed for in",0);
    util_file_unlink(in);
}

