      subroutine util_flush(unit)
C$Id: util_flush.F,v 1.22 2005-05-27 01:51:05 edo Exp $
      implicit none
#include "stdio.fh"
      integer unit
#ifdef CRAY
      integer*4 iunit4
#endif
c
c     Flush buffered output for the fortran unit
c
#if defined(SUN)      || defined(KSR)    || defined(PARAGON) || \
    defined(SOLARIS)  || \
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
#if defined(SGI)
      if (unit.eq.luout) then
        call sgi_flush6()
      else
        call flush(unit)
      endif
#endif
#if defined(IPSC) || defined(DELTA)
      call forflush(unit)
#endif
c
      end
