      subroutine util_flush(unit)
C$Id: util_flush.F,v 1.11 1995-06-21 18:29:02 pg511 Exp $
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
      call fflush(unit)      
#endif
#endif
#endif
c
      end
