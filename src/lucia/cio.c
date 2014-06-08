/* cio.c $Revision: 7.7 $ */
/******************************************************************************/
/*                                                                            */
/*                               A I X - I / O                                */
/*                                                                            */
/*  The fast I/O calls the following C-langue primitives:                     */
/*  open, close, read, write, lseek, remove and fsync                         */
/*  This file includes the FORTRAN to C-language interfaces.                  */

#include <fcntl.h>
#ifndef _WIN32_
#include <unistd.h>
#else
#include <windows.h>
#define open  _lopen
#define close _lclose
#define read  _lread
#define write _lwrite
#define lseek _llseek
#endif
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include "molcastype.h"

#define MIN(x,y) (x<y? x : y)
/*--- c_open -----------------------------------------------------------------*/

#ifdef _CAPITALS_
#define c_open C_OPEN
#else
#ifndef ADD_
#define c_open c_open_
#endif
#endif

INT c_open(Path)
 char *Path;

{

 INT rc;
 INT oFlag;
 INT oMode;
#ifdef _CRAY_C90_
 char fn[256];
 oFlag=O_CREAT|O_RDWR;
 (void)strcpy(fn,Path);
 rc = open(fn,oFlag,0644);
#else
#ifndef _WIN32_
 oFlag=O_CREAT|O_RDWR;
 oMode=S_IRUSR|S_IRGRP|S_IROTH|S_IWUSR;
 rc = open(Path,oFlag,oMode);
#else
 oFlag=OF_READWRITE;
 rc=open(Path,oFlag);
#endif
#endif
 if(rc<0) {
   oFlag=O_RDONLY;
#ifdef _CRAY_C90_
   rc = open(fn,oFlag);
#else
#ifndef _WIN32_
   rc = open(Path,oFlag);
#else
   oFlag=OF_READ;
   rc=open(Path,oFlag);
#endif
#endif
 }

 return rc;

}
/*--- c_open_w -----------------------------------------------------------------*/

#ifdef _CAPITALS_
#define c_openw C_OPENW
#else
#ifndef ADD_
#define c_openw c_openw_
#endif
#endif

INT c_openw(Path)
 char *Path;

{

 INT rc;
 INT oFlag;
 INT oMode;
#ifdef _CRAY_C90_
 char fn[256];
 oFlag=O_CREAT|O_RDWR|O_TRUNC;
 (void)strcpy(fn,Path);
 rc = open(fn,oFlag,0644);
#else
#ifndef _WIN32_
 oFlag=O_CREAT|O_RDWR|O_TRUNC;
 oMode=S_IRUSR|S_IRGRP|S_IROTH|S_IWUSR;
 rc = open(Path,oFlag,oMode);
#else
 oFlag=OF_READWRITE;
 rc=open(Path,oFlag);
#endif
#endif
 return rc;

}

/*--- c_close ----------------------------------------------------------------*/
#ifdef _CAPITALS_
#define c_close C_CLOSE
#else
#ifndef ADD_
#define c_close c_close_
#endif
#endif

INT c_close(FileDescriptor)
 INT *FileDescriptor;

{
 INT rc;
 rc = close(*FileDescriptor);
 return rc;
}

/*--- c_read -----------------------------------------------------------------*/

#ifdef _CAPITALS_
#define c_read C_READ
#else
#ifndef ADD_
#define c_read c_read_
#endif
#endif

INT c_read(FileDescriptor,Buffer,nBytes)
 INT *FileDescriptor;
 char *Buffer;
 INT *nBytes;

{
 INT rc=0;
 INT bfrblk=1024*1024;
 INT i=0;
 INT j;
 INT remains;
 INT readlength;
 remains=*nBytes;
 while (remains > 0){
      readlength = MIN(bfrblk,remains);
      rc = (INT)read(*FileDescriptor,(void *)(Buffer+i),(size_t)(readlength));
      if ( rc == readlength ) { i = i+readlength; rc = i; remains = remains - bfrblk;}
      else { rc = 0; return rc ;}
 }
 return rc;
}

/*--- c_write ----------------------------------------------------------------*/

#ifdef _CAPITALS_
#define c_write C_WRITE
#else
#ifndef ADD_
#define c_write c_write_
#endif
#endif

INT c_write(FileDescriptor,Buffer,nBytes)
 INT *FileDescriptor;
 char *Buffer;
 INT *nBytes;

{
 INT rc=0;
 INT bfrblk=1024*1024;
 INT i=0;
 INT remains;
 INT writelength;
 remains=*nBytes;
 while (remains > 0){
      writelength = MIN(bfrblk,remains);
      rc = (INT)write(*FileDescriptor,(void *)(Buffer+i),(size_t)(writelength));
      if ( rc == writelength ) { i = i+writelength; rc = i; remains = remains - bfrblk;}
      else { rc = 0; return rc ;}
 }
 return rc;
}

/*--- c_lseek ----------------------------------------------------------------*/

#ifdef _CAPITALS_
#define c_lseek C_LSEEK
#else
#ifndef ADD_
#define c_lseek c_lseek_
#endif
#endif

INT c_lseek(FileDescriptor,Offset)
 INT *FileDescriptor;
 INT *Offset;

{
#ifdef _WIN32_
typedef long off_t;
#endif
 INT rc;
 rc = (INT)lseek(*FileDescriptor,(off_t)(*Offset),SEEK_SET);
 return rc;
}

/*--- c_remove ---------------------------------------------------------------*/

#ifdef _CAPITALS_
#define c_remove C_REMOVE
#else
#ifndef ADD_
#define c_remove c_remove_
#endif
#endif

INT c_remove(FileName)
 char *FileName;

{
 INT rc;
#ifdef _CAPITALS_
 char fn[256];
#endif

#ifdef _CAPITALS_
 (void)strcpy(fn,FileName);
 rc = remove(fn);
#else
#ifndef _WIN32_
 rc = remove(FileName);
#else
 rc = DeleteFile(FileName);
#endif
#endif
 return rc;
}

/*--- c_fsync ----------------------------------------------------------------*/

#ifdef _CAPITALS_
#define c_fsync C_FSYNC
#else
#ifndef ADD_
#define c_fsync c_fsync_
#endif
#endif

INT c_fsync(FileDescriptor)
 INT *FileDescriptor;

{
 INT rc;
#ifndef _WIN32_
 rc = fsync(*FileDescriptor);
#else
 rc=0;
#endif
 return rc;
}

/*--- c_copy ----------------------------------------------------------------*/

#ifdef _CAPITALS_
#define c_copy C_COPY
#else
#ifndef ADD_
#define c_copy c_copy_
#endif
#endif

INT c_copy(FileDescriptor1, FileDescriptor2)
 INT *FileDescriptor1, *FileDescriptor2;
{
 INT rc;
 char *Buffer;
 struct stat stat;
 size_t rce;

 rc=fstat(*FileDescriptor1, &stat);

 rce=stat.st_size;
 Buffer=(char*) malloc(sizeof(char)*(rce+1));
      rc = (INT)read(*FileDescriptor1,(void *)(Buffer),(size_t)(rce));
      rc = (INT)write(*FileDescriptor2,(void *)(Buffer),(size_t)(rce));
 free(Buffer);
 return rc;
}
/* $Id$ */
