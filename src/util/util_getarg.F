      subroutine util_getarg(iarg,arg)
* $Id: util_getarg.F,v 1.2 1998-05-22 18:25:10 d3e129 Exp $
*
* The statement call util_getarg( k , arg ) gets the  kth  command-
* line argument and puts it into the string arg.
*
      implicit none
#if defined(SGI) || defined(SGITFP)
      integer*4 iarg
#else
      integer iarg
#endif
      character*(*) arg
#if defined(CRAY)
      integer length,ierror
#endif
*
#if defined(CRAY)
      call pxfgetarg(iarg,arg,length,ierror)
c
      if (ierror.ne.0)
     + call errquit('util_getarg:pxfgetarg fatal error',ierror)
#else
      call getarg(iarg,arg)
#endif
*
      end
