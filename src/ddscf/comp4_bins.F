      subroutine comp4_bins(m, n, nb_per_i)
C$Id$
      implicit none
c
      integer m(*)              ! Array holding info compressed to 4 bits
      integer n                 ! Bin the first n entries
      integer nb_per_i          ! No. of bytes per integer
c
#if defined(CRAY)
      integer vv,vvv            ! Temporaries for CRAY
#endif
      integer v                 ! Value after compression
      integer index, nbits, i
      integer bins(0:15)
      double precision fast(0:15)
#include "bitops_decls.fh"
      data fast/0.0d0, 1.0d-13, 1.0d-12, 1.0d-11, 1.0d-10, 1.0d-9,
     $     1.0d-8, 1.0d-7, 1.0d-6, 1.0d-5, 1.0d-4, 1.0d-3, 1.0d-2,
     $     1.0d-1, 1.0d0, 1.0d1/
#include "bitops_funcs.fh"
c
      do i = 0, 15
         bins(i) = 0
      enddo
c
      do i = 1, n
         v = 15
         index = (i - 1)/(2*nb_per_i) + 1
         nbits = 4*(i - (index-1)*(2*nb_per_i) - 1)
#if defined(CRAY)
         vvv =  ishft(v, nbits)
         vv = ishft(iand(m(index), vvv),-nbits)
         v = iand(vv,15)
#else
         v = iand(ishft(iand(m(index), ishft(v, nbits)), -nbits),15)
#endif
         bins(v) = bins(v) + 1
      enddo
c
      write(6,1) 
 1    format(1x,4('  value      no.   '))
      write(6,2) 
 2    format(1x,4(' -------- ------   '))
      write(6,3) (fast(i), bins(i), i=0,15)
 3    format(4(1x,4(1pd9.2,i7,3x)/))
c
      end
