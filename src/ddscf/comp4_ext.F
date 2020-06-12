      subroutine comp4_extract(m, i, s, nb_per_i)
C$Id$
      implicit none
c
      integer m(*)              ! Array holding info compressed to 4 bits
      integer i                 ! Index into compressed array
      double precision s        ! Returns extracted value
      integer nb_per_i          ! No. of bytes per integer
c
      integer v                 ! Value after compression
#if defined(CRAY)
      integer vv, vvv
#endif
      integer index, nbits
#include "bitops_decls.fh"
      double precision fast(0:15)
      data fast/0.0d0, 1.0d-13, 1.0d-12, 1.0d-11, 1.0d-10, 1.0d-9,
     $     1.0d-8, 1.0d-7, 1.0d-6, 1.0d-5, 1.0d-4, 1.0d-3, 1.0d-2,
     $     1.0d-1, 1.0d0, 1.0d1/
#include "bitops_funcs.fh"
c
      v = 15
      index = (i - 1)/(2*nb_per_i) + 1
      nbits = 4*(i - (index-1)*(2*nb_per_i) - 1)
#if defined(CRAY)
      vvv =  shiftl(v, nbits)
      vv = shiftr(iand(m(index), vvv), nbits)
      v = iand(vv,15)
#else
      v = iand(ishft(iand(m(index), ishft(v, nbits)), -nbits),15)
#endif
c
      s = fast(v)
*      write(6,*) i, ' -> ', index, nbits, v, s
c
      end
