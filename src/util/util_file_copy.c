/*
 $Id: util_file_copy.c,v 1.12 2003-08-13 18:06:11 edo Exp $
 */

#include <stdio.h>
#if defined(CRAY) && !defined(__crayx1)
#include <fortran.h>
#define FATR
#endif
#include "global.h"
#include "typesf2c.h"

#if defined(USE_FCD)
int fortchar_to_string(_fcd, int, char *, const int);
#else
Integer fortchar_to_string(const char *, Integer, char *, const Integer);
#endif

void util_file_copy(const char *input, const char *output)
/*
  The local process copies the input file to the output file.
  Any existing file is overwritten.
  Any error is fatal.
  */
{
    FILE *fin  = fopen(input, "rb");
    FILE *fout = fopen(output, "w+b");
    char buf[8192];
    Integer nread;

    if (!fin) {
	fprintf(stderr,"util_file_copy: unable to open %s\n", input);
	ga_error("util_file_copy",0);
    }
    if (!fout) {
	fprintf(stderr,"util_file_copy: unable to open %s\n", output);
	ga_error("util_file_copy",0);
    }
    while ((nread = fread(buf, 1, sizeof(buf), fin)) > 0)
	if (fwrite(buf, 1, nread, fout) != nread) {
	    fprintf(stderr,"util_file_copy: failed writing %s\n", output);
	    ga_error("util_file_copy",0);
	}

    if (!feof(fin)) {
	fprintf(stderr,"util_file_copy: failed reading %s\n", input);
	ga_error("util_file_copy",0);
    }
	
    (void) fclose(fin);
    (void) fclose(fout);
}

void util_file_parallel_copy(const char *input, const char *output)
/*
  Process 0 broadcasts the input file to all nodes which write to their
  output file.  Any existing file is overwritten.  Any error is fatal.

  It is same to call it with a same/different output file name on all nodes.
  If the input and output files have the same name on node0 it is also safe.
  */
{
    FILE *fin=0, *fout=0;
    Integer differ = strcmp(input,output);

    if (ga_nodeid_() == 0) {
      if (!(fin = fopen(input, "rb"))) {
	fprintf(stderr,"util_file_copy: unable to open input %s\n", input);
	ga_error("util_file_parallel_copy",0);
      }
      if (differ) {
	if (!(fout = fopen(output, "w+b"))) {
	  fprintf(stderr,"util_file_copy: unable to open output %s\n", input);
	  ga_error("util_file_parallel_copy",0);
	}
      }
    }
    else if (!(fout = fopen(output, "w+b"))) {
      fprintf(stderr,"util_file_copy: unable to open output %s\n", input);
      ga_error("util_file_parallel_copy",0);
    }

    while (1) {
      char buf[8192];
      Integer nread, msgnread=44,msgbuf=45,msglen=sizeof(Integer),node0=0;
      if (ga_nodeid_() == 0)
	nread = fread(buf, 1, sizeof(buf), fin);
      ga_brdcst_(&msgnread, &nread, &msglen, &node0);
      if (nread > 0) {
	ga_brdcst_(&msgbuf, buf, &nread, &node0);
	if ((ga_nodeid_() != 0) || (differ != 0)) {
	  if (fwrite(buf, 1, nread, fout) != nread) {
	    fprintf(stderr,"util_file_parallel_copy: failed writing %s\n", output);
	    ga_error("util_file_parallel_copy",0);
	  }
	}
      }
      else
	break;
    }
	
    if (ga_nodeid_() == 0) {
      if (!feof(fin)) {
	fprintf(stderr,"util_file_parallel_copy: failed reading %s\n", input);
	ga_error("util_file_parallel_copy",0);
      }
    }
	
    if (fin) (void) fclose(fin);
    if (fout) (void) fclose(fout);
}

#if defined(USE_FCD)
void FATR UTIL_FILE_COPY(_fcd input, _fcd output)
{
    int lin  = _fcdlen(input);
    int lout = _fcdlen(output);
#else
void util_file_copy_(const char *input, const char *output, Integer lin, Integer lout)
{
#endif
    char in[255], out[255];
    if (!fortchar_to_string(input, lin, in, sizeof(in)))
	ga_error("util_file_copy: fortchar_to_string failed for in",0);
    if (!fortchar_to_string(output, lout, out, sizeof(out)))
	ga_error("util_file_copy: fortchar_to_string failed for out",0);
    util_file_copy(in, out);
}

