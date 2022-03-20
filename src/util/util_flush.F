      subroutine util_flush(unit)
C$Id$
      implicit none
#include "stdio.fh"
      integer unit
#ifdef CRAY
      integer*4 iunit4
#endif
c
c     Flush buffered output for the fortran unit
c
#if defined(SUN) || defined(SOLARIS) || \
    defined(IBM)||defined(LINUX)||defined(MACX)
      call flush(unit)
#endif
#ifdef CRAY
      iunit4=unit
      call flush(iunit4)
#endif
#if defined(PSCALE)
      call util_cflush()
#endif
#if defined(GFORTRAN)
      call flush(unit)
#endif
      call util_cflush()
c
      end
