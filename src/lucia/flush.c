/* RCS *******************************************************************/
/* $Id$ */
/* $Log: flush.c,v $
/* Revision 1.1  2004/06/16 19:27:22  andreas
/* new routines and includes for
/*  -- densities with unrestricted orbitals
/*  -- a file handler
/*  -- a non-linear optimization controller
/*  -- orbital-optimized CC
/*
/* Revision 1.4  2003/10/13 12:11:33  dimi
/* Bug on SGI/Irix fixed
/*
/* Revision 1.3  2003/01/21 22:14:39  uwe
/* flush under Linux too
/*
/* Revision 1.2  2000/10/20 13:19:39  haettig
/* removed some `/' in the  comment lines of the file headers.
/*
 * Revision 1.1  1992/09/17 10:56:23  marco
 * Initial revision
 * */
/* RCS *******************************************************************/

#include <stdio.h>

#if !defined(IRIS) && !defined(CONVEX)
int flush(iunit)
   int *iunit;
{
   if (*iunit == 6) fflush(stdout);
}
#endif

/* Intel Compiler 6.0 and 7.0 have problems wiht flushing output if it is appended to an existing file..., so lets do it in C*/
#if defined(LINUX) || defined(IRIS)
int flush_(iunit)
   int *iunit;
{
   if (*iunit == 6) fflush(stdout);
}
#endif
