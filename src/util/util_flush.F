      subroutine util_flush(unit)
C$Id: util_flush.F,v 1.9 1995-02-07 00:15:32 rg240 Exp $
      implicit none
      integer unit
c
c     Flush buffered output for the fortran unit
c
#if defined(SUN) || defined(KSR)    || defined(PARAGON) || \
    defined(SGI) || defined(SGITFP) || defined(IBM) || \
    defined(CRAY_T3D)
      call flush(unit)
#elif defined(IPSC) || defined(DELTA)
      call forflush(unit)
#endif
c
      end
