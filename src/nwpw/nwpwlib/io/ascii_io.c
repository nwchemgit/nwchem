/*
 $Id: ascii_io.c 26153 2014-09-03 05:21:19Z edo $
*********************************************************************************
*									        *
* ascii_io.c		   					                *
*									        *
* These routines implement a simple input/output package for writing to	        *
* and reading from gzipped files.					        *
*                                                                               *
*									        *
* Fortran usage:							        *
*	call ascii_openfile(unit,'file', 'X', n)	X is either 'r' or 'w'	*
*					n is the length of the filename	        *
*	call ascii_cwrite(unit,c, n)		c is n element character array	*
*	call ascii_iwrite(unit,i, n)		i is n element integer array	*
*	call ascii_dwrite(unit,d, n)		d is n element double array	*
*	call ascii_cread(unit,c, n)		c is n element character array	*
*	call ascii_iread(unit,i, n)		i is n element integer array	*
*	call ascii_dread(unit,d, n)		d is n element double array	*
*	call ascii_closefile(unit)		close the datafile		*
*									        *
* Author:      Scott Kohn (skohn@chem.ucsd.edu)				        *
* modified by: Eric Bylaska (ebylaska@chem.ucsd.edu)			        *
*									        *
*********************************************************************************
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "typesf2c.h"

#if !defined(__crayx1)

#if (defined(CRAY) || defined(WIN32)) && !defined(__crayx1) &&!defined(__MINGW32__)
#define ascii_cwrite_ ASCII_CWRITE
#define ascii_cread_  ASCII_CREAD
#define ascii_iwrite_ ASCII_IWRITE
#define ascii_iread_  ASCII_IREAD
#define ascii_dwrite_ ASCII_DWRITE
#define ascii_dread_  ASCII_DREAD
#define ascii_dshift_fileptr_  ASCII_DSHIFT_FILEPTR
#define ascii_ishift_fileptr_  ASCII_ISHIFT_FILEPTR
#define ascii_openfile_  ASCII_OPENFILE
#define ascii_closefile_ ASCII_CLOSEFILE
#endif

#endif


#define MAX_UNIT	10

static FILE* fd[MAX_UNIT];	/* the file descriptor of the pipe */

#define BAIL(X) { fprintf(stderr, X); exit(-1); }

/*
*************************************************************************
*									*
* Define the Xwrite and Xread routines using the Fortran bindings.	*
*									*
*************************************************************************
*/

void FATR ascii_cwrite_
#if defined(USE_FCD)
(const Integer *unit, _fcd fcd_c, const Integer *n)
{
   const char *c = _fcdtocp(fcd_c);

#else
(const Integer *unit, char *c, const Integer *n)
{
#endif

   int j; 
   for (j=0; j<(*n); ++j)
      (void) fprintf(fd[*unit],"%c ",c[j]);
   (void) fprintf(fd[*unit],"\n");
}

void FATR ascii_cread_
#if defined(USE_FCD)
(const Integer *unit, _fcd fcd_c, const Integer *n)
{
    char *c = _fcdtocp(fcd_c);

#else
(const Integer *unit, char *c, const Integer *n)
{
#endif

   int j; 
   for (j=0; j<(*n); ++j)
      (void) fscanf(fd[*unit],"%c ",&c[j]);
}

void FATR ascii_iwrite_(const Integer *unit, const Integer *i, const Integer *n)
{
   int j;
   for (j=0; j<(*n); ++j)
      fprintf(fd[*unit],"%d ",(int) i[j]);
   fprintf(fd[*unit],"\n");

}


void FATR ascii_iread_(const Integer *unit, Integer *i, const Integer *n)
{
   int j,k;
   for (j=0; j<(*n); ++j)
   {
      (void) fscanf(fd[*unit],"%d ",&k);
      i[j] = (Integer) k;
   }
}

void FATR ascii_dwrite_(const Integer *unit, const DoublePrecision *d, const Integer *n)
{
   int j;
   for (j=0; j<(*n); ++j)
      fprintf(fd[*unit],"%20.15le ",d[j]);
   fprintf(fd[*unit],"\n");

}

void FATR ascii_dread_(const Integer *unit, DoublePrecision *d, const Integer *n)
{
   int j;
   for (j=0; j<(*n); ++j)
      (void) fscanf(fd[*unit],"%le ",&d[j]);
}

void FATR ascii_dshift_fileptr_(const Integer *unit, const Integer *n)
{
   (void) fseek(fd[*unit], ((long) (*n)*sizeof(DoublePrecision)),SEEK_CUR);
}
void FATR ascii_ishift_fileptr_(const Integer *unit, const Integer *n)
{
   (void) fseek(fd[*unit], ((long) (*n)*sizeof(Integer)),SEEK_CUR);
}

/*
*************************************************************************
*									*
* void ascii_openfile(char *filename, char *mode, Integer *n)		*
* void ascii_closefile()						*
*									*
* Function openfile opens a pipe to either gzip (to compress a stream)	*
* or zcat (to uncompress a stream).  Function ascii_closefile() closes 	*
* the pipe stream created by ascii_openfile().				*
*									*
*************************************************************************
*/

#define FUDGE_FACTOR (8)

void FATR ascii_openfile_
#if defined(USE_FCD)
(const Integer *unit, _fcd fcd_filename, Integer *n1, _fcd fcd_mode, Integer *n2)
{
   const char *filename = _fcdtocp(fcd_filename);
   const char *mode = _fcdtocp(fcd_mode);

#else
(const  Integer *unit, char *filename, Integer *n1, char *mode,    Integer *n2)
{
#endif

   char *file = (char *) malloc(*n1+1);

   (void) strncpy(file, filename, *n1);
   file[*n1] = '\0';

   if ((*mode == 'r') || (*mode == 'R')) {
      if (!(fd[*unit] = fopen(file, "r")))
         BAIL("ERROR:  Could not open pipe from input file\n");
   } else {
      if (!(fd[*unit] = fopen(file, "w")))
         BAIL("ERROR:  Could not open pipe to output file\n");
   }

   free(file);

}

void FATR ascii_closefile_(const Integer *unit)
{
   (void) fclose(fd[*unit]);
}


