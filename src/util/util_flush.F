      subroutine util_flush(unit)
C$Id: util_flush.F,v 1.13 1996-01-10 15:58:46 d3g681 Exp $
      implicit none
      integer unit
c
c     Flush buffered output for the fortran unit
c
#if defined(SUN) || defined(KSR)    || defined(PARAGON) || \
    defined(SGI) || defined(SGITFP) || \
    defined(CRAY_T3D) || defined(IBM)
      call flush(unit)
#endif
#if defined(IPSC) || defined(DELTA)
      call forflush(unit)
#endif
c
      end
