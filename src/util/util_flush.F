      subroutine util_flush(unit)
C$Id: util_flush.F,v 1.15 1997-03-04 06:07:10 d3e129 Exp $
      implicit none
      integer unit
c
c     Flush buffered output for the fortran unit
c
#if defined(SUN)      || defined(KSR)    || defined(PARAGON) || \
    defined(SGI)      || defined(SGITFP) || defined(SOLARIS) || \
    defined(CRAY_T3D) || defined(CRAY_T3E) || defined(IBM)
      call flush(unit)
#endif
#if defined(IPSC) || defined(DELTA)
      call forflush(unit)
#endif
c
      end
