      subroutine comp4_insert(m, i, s, nb_per_i)
C$Id$
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
c     Map value to log10 with range 1,...,-14 and values outside
c     of this range mapped to nearest limit.
c
c     Previously rounded values inside range to nearest integer log value,
c     however, since we use the data for screening, it's better to
c     always round up so that we preserve the upper bound properties
c     of the screening.  BUT, this leads to a possible overestimate of
c     (ij|kl) by a factor of 100, whereas rounding to the nearest is
c     just off by a factor of 10 either way. For very sparse cases
c     this can be significant so stick with the round to nearest.
c
      if (abs(s) .le. 1.0d-14) then
         v = -14
      else if (abs(s) .ge. 1.0d1) then
         v = 1
      else
      v = nint(log10(abs(s)))  ! round to nearest
*     v = int(log10(abs(s)) + 0.99d0) ! round up
      endif
c
c     Finally shift +1,...,-14, to 15,...,0
c
      v = v + 14
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
      m(index) = ior(m(index), ishft(v, nbits))
#endif
c
      end
