      subroutine util_flush(unit)
C$Id: util_flush.F,v 1.8 1995-02-02 18:09:53 d3g681 Exp $
      implicit none
      integer unit
c
c     Flush buffered output for the fortran unit
c
#if defined(SUN) || defined(KSR)    || defined(PARAGON) || \
    defined(SGI) || defined(SGITFP) || defined(IBM)
      call flush(unit)
#elif defined(IPSC) || defined(DELTA)
      call forflush(unit)
#endif
c
      end
