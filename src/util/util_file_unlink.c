/*
 $Id: util_file_unlink.c,v 1.3 1999-10-21 22:54:44 d3g681 Exp $
 */

/*
  Try to route all sequential file operations thru eaf/elio
  so that can handle extents correctly
*/


#include <stdio.h>
#include <unistd.h>
#ifdef CRAY
#include <fortran.h>
#endif

#include "eaf.h"

#ifdef CRAY
int fortchar_to_string(_fcd, int, char *, const int);
#else
int fortchar_to_string(const char *, int, char *, const int);
#endif

void ga_error(const char *, long);

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
	    ga_error("util_file_unlink",0);
	}
    }
  */
  if (eaf_delete(filename) != 0)
    ga_error("util_file_unlink",0);
}

#ifdef CRAY
void UTIL_FILE_UNLINK(_fcd input)
{
    int lin  = _fcdlen(input);
#else
void util_file_unlink_(const char *input, int lin)
{
#endif
    char in[255];
    if (!fortchar_to_string(input, lin, in, sizeof(in)))
	ga_error("util_file_unlink: fortchar_to_string failed for in",0);
    util_file_unlink(in);
}

