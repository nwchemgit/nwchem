      subroutine util_flush(unit)
C$Id: util_flush.F,v 1.12 1995-06-23 19:30:32 pg511 Exp $
      implicit none
      integer unit
c
c     Flush buffered output for the fortran unit
c
#if defined(SUN) || defined(KSR)    || defined(PARAGON) || \
    defined(SGI) || defined(SGITFP) || \
    defined(CRAY_T3D)
      call flush(unit)
#else
#if defined(IPSC) || defined(DELTA)
      call forflush(unit)
#else
#if defined(IBM)

      return
#endif
#endif
#endif
c
      end
