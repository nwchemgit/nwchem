      subroutine comp4_insert(m, i, s, nb_per_i)
C$Id: comp4_insert.F,v 1.6 1995-11-03 20:41:53 d3g681 Exp $
      implicit none
c
      integer m(*)              ! Array holding info compressed to 4 bits
      integer i                 ! Index into compressed array
      double precision s        ! Value to be inserted
      integer nb_per_i          ! No. of bytes per integer
c
      integer v                 ! Value after compression
      integer index, nbits      ! Index into compressed array
#ifdef CRAY
      integer vv
#endif
#include "bitops.fh"      
c
c     Map value to log10 with range 4,...,-11 and values outside
c     of this range mapped to nearest limit
c
      if (abs(s) .le. 1.0d-11) then
         v = -11
      else if (abs(s) .ge. 1.0d4) then
         v = 4
      else
         v = nint(log10(abs(s)))
      endif
c
c     Finally shift +4,...,-11, to 15,...,0
c
      v = v + 11
c
c     Now which integer holds this and how many bits must v be shifted
c     before ORing in the value
c
      index = (i - 1)/(2*nb_per_i) + 1
      nbits = 4*(i - (index-1)*(2*nb_per_i) - 1)
*      write(6,*) s, ' -> ', v, index, nbits
c
#if defined(CRAY)
      vv=shiftl(v, nbits)
      m(index) = ior(m(index),vv)
#else
      m(index) = ior(m(index), lshift(v, nbits))
#endif
c
      end
