      subroutine util_flush(unit)
      implicit none
      integer unit
c
c     Flush buffered output for the fortran unit
c
#if defined(SUN) || defined(KSR)
      call flush(unit)
#elif defined(IPSC)
      call forflush(unit)
#endif
c
      end
