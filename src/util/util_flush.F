      subroutine util_flush(unit)
C$Id: util_flush.F,v 1.16 1999-07-02 00:30:12 d3e129 Exp $
      implicit none
#include "stdio.fh"
      integer unit
c
c     Flush buffered output for the fortran unit
c
#if defined(SUN)      || defined(KSR)    || defined(PARAGON) || \
    defined(SOLARIS)  || \
    defined(CRAY_T3D) || defined(CRAY_T3E) || defined(IBM)
      call flush(unit)
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
