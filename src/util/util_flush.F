      subroutine util_flush(unit)
      implicit none
      integer unit
c
c     Flush buffered output for the fortran unit
c
#if defined(SUN) || defined(KSR) || defined(PARAGON) || \
    defined(SGI) || defined(SGITFP)
      call flush(unit)
#elif defined(IPSC) || defined(DELTA)
      call forflush(unit)
#endif
c
      end
