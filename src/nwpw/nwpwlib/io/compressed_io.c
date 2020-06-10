/*
 $Id$
*************************************************************************
*									*
* compressed_io.c							*
*									*
* These routines implement a simple input/output package for writing to	*
* and reading from gzipped files.					*
* If filename ends in .ascii then ascii io is used                      *
*									*
* Fortran usage:							*
*	call openfile(unit,'file', 'X', n)	X is either 'r' or 'w' 	*
*					n is the length of the filename	*
*	call cwrite(unit,c, n)		c is n element character array	*
*	call iwrite(unit,i, n)		i is n element integer array	*
*	call dwrite(unit,d, n)		d is n element double array	*
*	call cread(unit,c, n)		c is n element character array	*
*	call iread(unit,i, n)		i is n element integer array	*
*	call dread(unit,d, n)		d is n element double array	*
*	call closefile(unit)		close the datafile		*
*									*
* Author:      Scott Kohn (skohn@chem.ucsd.edu)				*
* modified by: Eric Bylaska (ebylaska@chem.ucsd.edu)			*
*									*
*************************************************************************
*/

#include <stdio.h>
#include <unistd.h>

#include <stdlib.h>
#include <string.h>
#include "typesf2c.h"

#if !defined(__crayx1)

#if defined(CRAY) || defined(CRAY_T3D) 
#include <fortran.h>
#define USE_FCD
#endif

#if (defined(CRAY) || defined(WIN32)) && !defined(__crayx1) &&!defined(__MINGW32__)
#define cwrite_ CWRITE
#define cread_  CREAD
#define iwrite_ IWRITE
#define iread_  IREAD
#define dwrite_ DWRITE
#define dread_  DREAD
#define dshift_fileptr_  DSHIFT_FILEPTR
#define ishift_fileptr_  ISHIFT_FILEPTR
#define flush_fileptr_  FLUSH_FILEPTR
#define openfile_  OPENFILE
#define closefile_ CLOSEFILE
#endif

#endif


#define MAX_UNIT	10

static FILE* fd[MAX_UNIT];	/* the file descriptor of the pipe */

static int ascii[MAX_UNIT];

#define BAIL(X) { fprintf(stderr, X); exit(-1); }

/*
*************************************************************************
*									*
* Define the Xwrite and Xread routines using the Fortran bindings.	*
*									*
*************************************************************************
*/

void FATR cwrite_
#if defined(USE_FCD)
(const Integer *unit, 
                   _fcd fcd_c,
             const Integer *n)
{
   const char *c = _fcdtocp(fcd_c);

#else
(const Integer *unit, 
                                char *c,
                          const Integer *n)
{
#endif

   if (ascii[*unit])
   {
      int j;
      for (j=0; j<(*n); ++j)
         (void) fprintf(fd[*unit],"%c ",c[j]);
      (void) fprintf(fd[*unit],"\n");
   }
   else
      (void) fwrite(c, sizeof(char), *n, fd[*unit]);
}

void FATR cread_
#if defined(USE_FCD)
(const Integer *unit, 
                  _fcd fcd_c,
            const Integer *n)
{
    char *c = _fcdtocp(fcd_c);

#else
(const Integer *unit, 
                     char *c, 
            const Integer *n)
{
#endif

   if (ascii[*unit])
   {
      int j;
      for (j=0; j<(*n); ++j)
         (void) fscanf(fd[*unit],"%c ",&c[j]);
   }
   else
      (void) fread(c, sizeof(char), *n, fd[*unit]);
}

void FATR iwrite_(const Integer *unit, const Integer *i, const Integer *n)
{
   if (ascii[*unit])
   {
      int j;
      for (j=0; j<(*n); ++j)
         fprintf(fd[*unit],"%d ",(int) i[j]);
      fprintf(fd[*unit],"\n");
   }
   else
      (void) fwrite(i, sizeof(Integer), *n, fd[*unit]);
}

void FATR iread_(const Integer *unit, Integer *i, const Integer *n)
{
   if (ascii[*unit])
   {
      int j,k;
       for (j=0; j<(*n); ++j)
       {
          (void) fscanf(fd[*unit],"%d ",&k);
          i[j] = (Integer) k;
       }
   }
   else
      (void) fread(i, sizeof(Integer), *n, fd[*unit]);
}

void FATR dwrite_(const Integer *unit, const DoublePrecision *d, const Integer *n)
{
   if (ascii[*unit])
   {
      int j;
      for (j=0; j<(*n); ++j)
         fprintf(fd[*unit],"%20.15le ",d[j]);
      fprintf(fd[*unit],"\n");
   }
   else
      (void) fwrite(d, sizeof(DoublePrecision), *n, fd[*unit]);
}

void FATR dread_(const Integer *unit, DoublePrecision *d, const Integer *n)
{
   if (ascii[*unit])
   {
      int j;
      for (j=0; j<(*n); ++j)
         (void) fscanf(fd[*unit],"%le ",&d[j]);
   }
   else
      (void) fread(d, sizeof(DoublePrecision), *n, fd[*unit]);
}

void FATR dshift_fileptr_(const Integer *unit, const Integer *n)
{
   (void) fseek(fd[*unit], ((long) (*n)*sizeof(DoublePrecision)),SEEK_CUR);
}
void FATR ishift_fileptr_(const Integer *unit, const Integer *n)
{
   (void) fseek(fd[*unit], ((long) (*n)*sizeof(Integer)),SEEK_CUR);
}
void FATR flush_fileptr_(const Integer *unit)
{
   (void) fflush(fd[*unit]);
#ifdef WIN32   
   (void) _commit(fileno(fd[*unit]));
#else
   (void) fsync(fileno(fd[*unit]));
#endif   
}

/*
*************************************************************************
*									*
* void openfile(char *filename, char *mode, Integer *n)			*
* void closefile()							*
*									*
* Function openfile opens a pipe to either gzip (to compress a stream)	*
* or zcat (to uncompress a stream).  Function closefile() closes the	*
* pipe stream created by openfile().					*
*									*
*************************************************************************
*/

#define FUDGE_FACTOR (8)

void FATR openfile_
#if defined(USE_FCD)
(const Integer *unit, 
                      _fcd fcd_filename,
                      Integer *n1,
       		      _fcd fcd_mode,
		      Integer *n2)
{
   const char *filename = _fcdtocp(fcd_filename);
   const char *mode = _fcdtocp(fcd_mode);

#else
(const  Integer *unit, 
                      char *filename, Integer *n1,
       		      char *mode,    
		      Integer *n2)
{
#endif

   char *file = (char *) malloc(*n1+1);
   char *rc, *wc;


   
   (void) strncpy(file, filename, *n1);
   file[*n1] = '\0';

   /* look for ".ascii" substring */
   if (strstr(filename,".ascii"))
   {
      ascii[*unit] = 1;
      rc = "r";
      wc = "w";
   }
   else
   {
      ascii[*unit] = 0;
      rc = "rb";
      wc = "wb";
   }

   if ((*mode == 'r') || (*mode == 'R')) {
      if (!(fd[*unit] = fopen(file, rc)))
         BAIL("ERROR:  Could not open pipe from input file\n");
   } 
   else
   {
      if (!(fd[*unit] = fopen(file, wc)))
         BAIL("ERROR:  Could not open pipe to output file\n");
   }

   free(file);

}

void FATR closefile_(const Integer *unit)
{
   (void) fclose(fd[*unit]);
}


