      subroutine util_getarg(iarg,arg)
* $Id: util_getarg.F,v 1.1 1998-05-18 21:28:06 d3j191 Exp $
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
      call getenv(iarg,arg)
#endif
*
      end
