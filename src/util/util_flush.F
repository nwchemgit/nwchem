      subroutine util_flush(unit)
C$Id: util_flush.F,v 1.10 1995-06-02 05:24:53 gg502 Exp $
      implicit none
      integer unit
c
c     Flush buffered output for the fortran unit
c
#if defined(SUN) || defined(KSR)    || defined(PARAGON) || \
    defined(SGI) || defined(SGITFP) || defined(IBM) || \
    defined(CRAY_T3D)
      call flush(unit)
#else
#if defined(IPSC) || defined(DELTA)
      call forflush(unit)
#endif
#endif
c
      end
