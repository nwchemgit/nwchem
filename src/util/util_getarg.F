      subroutine util_getarg(iarg_in,arg)
* $Id$
*
* The statement call util_getarg( k , arg ) gets the  kth  command-
* line argument and puts it into the string arg.
*
      implicit none
#include "errquit.fh"
      integer iarg_in
      character*(*) arg
#if defined(USE_FORTRAN2003)
      call get_command_argument(iarg_in,arg)
#else
#if defined(EXT_INT) && !( __GNUC__ >= 4)
      integer*4 iarg
#else
      integer iarg
#if defined(WIN32) &&!defined(__MINGW32__)
      integer idum
#endif
#endif
#if defined(CRAY)
      integer*4 length,ierror
#endif
#if defined(HPUX)
#ifdef EXT_INT
!$HP$ ALIAS getarg='getarg__'
#else
!$HP$ ALIAS getarg='getarg_'
#endif
#endif
*
      iarg=iarg_in
#if defined(CRAY)
      call pxfgetarg(iarg,arg,length,ierror)
c
      if (ierror.ne.0)
     + call errquit('util_getarg:pxfgetarg fatal error',ierror,
     &       UNKNOWN_ERR)
#elif defined(WIN32) &&!defined(__MINGW32__)
      call getarg(iarg,arg,idum)
#else
      call getarg(iarg,arg)
#endif
*
#endif
      return
      end
