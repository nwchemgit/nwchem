      subroutine util_getarg(iarg_in,arg)
* $Id: util_getarg.F,v 1.4 2000-08-14 23:04:51 bjohnson Exp $
*
* The statement call util_getarg( k , arg ) gets the  kth  command-
* line argument and puts it into the string arg.
*
      implicit none
      integer iarg_in
      character*(*) arg
#if defined(EXT_INT) && !defined(HPUX64)
      integer*4 iarg
#else
      integer iarg
#ifdef WIN32
cDEC$ IF (_DF_VERSION_ .lt. 600)
      integer idum
cDEC$ ENDIF
#endif
#endif
#if defined(CRAY)
      integer length,ierror
#endif
#if defined(HPUX64)
C
C    this hack allows to get the f77 external getarg out of libU77
C    insted of the f90 intrinsic
C
      external getarg
#endif
*
      iarg=iarg_in
#if defined(CRAY)
      call pxfgetarg(iarg,arg,length,ierror)
c
      if (ierror.ne.0)
     + call errquit('util_getarg:pxfgetarg fatal error',ierror)
#else
#ifdef WIN32
cDEC$ IF (_DF_VERSION_ .lt. 600)
      call getarg(iarg,arg,idum)
cDEC$ ELSE
      call getarg(iarg,arg)
cDEC$ ENDIF
#else
      call getarg(iarg,arg)
#endif
#endif
*
      return
      end
