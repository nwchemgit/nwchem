#include <stdio.h>
#ifdef CRAY
#include <fortran.h>
#endif

#ifdef CRAY
int fortchar_to_string(_fcd, int, char *, const int);
#else
int fortchar_to_string(const char *, int, char *, const int);
#endif

void ga_error(const char *, long);

void util_file_copy(const char *input, const char *output)
/*
  The local process copies the input file to the output file.
  Any existing file is overwritten.
  Any error is fatal.
  */
{
    FILE *fin  = fopen(input, "r");
    FILE *fout = fopen(output, "w+");
    char buf[8192];
    int nread;

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

#ifdef CRAY
void UTIL_FILE_COPY(_fcd input, _fcd output)
{
    int lin  = _fcdlen(input);
    int lout = _fcdlen(output);
#else
void util_file_copy_(const char *input, const char *output, int lin, int lout)
{
#endif
    char in[255], out[255];
    if (!fortchar_to_string(input, lin, in, sizeof(in)))
	ga_error("util_file_copy: fortchar_to_string failed for in",0);
    if (!fortchar_to_string(output, lout, out, sizeof(out)))
	ga_error("util_file_copy: fortchar_to_string failed for out",0);
    util_file_copy(in, out);
}

