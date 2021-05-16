      subroutine util_pack_info(bytes_per_int, obigendian)
      implicit none
#if defined(LINUX) || defined(LINUX64)
#define INTEGER_1 integer*1
#else
#define INTEGER_1 byte
#endif
#include "mafdecls.fh"
      integer bytes_per_int
      logical obigendian
c
      integer longw
      integer*2 shortw(4),mold
#ifdef F77      
      equivalence (longw,shortw)
#endif      
c     
      bytes_per_int = ma_sizeof(mt_int, 1, mt_byte)
c
      shortw(1) = 0
      longw = 1
#ifndef F77
      shortw=transfer(longw,mold,bytes_per_int/2)
#endif      
      obigendian = .not. (shortw(1) .eq. 1)
c
      end
#ifdef BYTE_SHORT_PACKING
      subroutine util_pack_8(nunpacked, packed, unpacked)
      implicit none
      integer nunpacked
      INTEGER_1 packed(*), unpacked(*)
c
      integer i, i1, bytes_per_int
      logical obigendian
c
      call util_pack_info(bytes_per_int, obigendian)
c
      if (obigendian) then
         i1 = bytes_per_int
      else
         i1 = 1
      endif
c
      do i = 1, nunpacked
         packed(i) = unpacked(i1)
         i1 = i1 + bytes_per_int
      enddo
c
      end
      subroutine util_unpack_8(nunpacked, packed, unpacked)
      implicit none
      integer nunpacked
      INTEGER_1 packed(*), unpacked(*)
c
      integer i, i1, ilo, ihi, ndo, base
      integer bytes_per_int
      logical obigendian
c
      call util_pack_info(bytes_per_int, obigendian)
c
      if (obigendian) then
         i1   = bytes_per_int
         base = bytes_per_int - 1
      else
         i1   = 1
         base = 0
      endif
c
c     This routine has been tweaked for IBM P2SC for which
c     the slow clock and high memory bandwidth does not
c     favour bit operations in the inner loop.
c
      do ilo = 1, nunpacked, 32768       ! 128 K cache
         ihi = min(nunpacked, ilo+32768-1)
         ndo = ihi - ilo + 1
         call ifill(ndo, 0, unpacked(i1-base), 1)
         do i = ilo, ihi
            unpacked(i1) = packed(i)
            i1 = i1 + bytes_per_int
         enddo
      enddo
c
      end
      subroutine util_pack_16(nunpacked, packed, unpacked)
*
* $Id$
*
      implicit none
      integer nunpacked
      integer*2 packed(*), unpacked(*)
c
      integer i, i2
      integer bytes_per_int, shorts_per_int
      logical obigendian
c
      call util_pack_info(bytes_per_int, obigendian)
      shorts_per_int = bytes_per_int/2
c
      if (obigendian) then
         i2 = shorts_per_int
      else
         i2 = 1
      endif
      do i = 1, nunpacked
         packed(i) = unpacked(i2)
         i2 = i2 + shorts_per_int
      enddo
c
      end
      subroutine util_unpack_16(nunpacked, packed, unpacked)
      implicit none
      integer nunpacked
      integer*2 packed(*), unpacked(*)
c
      integer i, i2, ilo, ihi, ndo, base
      integer bytes_per_int, shorts_per_int
      logical obigendian
c
      call util_pack_info(bytes_per_int, obigendian)
      shorts_per_int = bytes_per_int/2
c
      if (obigendian) then
         i2   = shorts_per_int
         base = shorts_per_int - 1
      else
         i2   = 1
         base = 0
      endif
c
      do ilo = 1, nunpacked, 32768       ! 128 K cache
         ihi = min(nunpacked, ilo+32768-1)
         ndo = ihi - ilo + 1
         call ifill(ndo, 0, unpacked(i2-base), 1)
         do i = ilo, ihi
            unpacked(i2) = packed(i)
            i2 = i2 + shorts_per_int
         enddo
      enddo
c
      end
#else
      subroutine util_pack_8(nunpacked, packed, unpacked)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
c
c     pack 8 bit integers ... n must be a multiple of 8
c
      integer nunpacked, packed(*), unpacked(*)
      integer npacked_per_int
      integer i, ij, kl, ijkl, ind, mn, pq, mnpq
#include "bitops.fh"
c
      npacked_per_int = ma_sizeof(mt_int, 1, mt_byte)
c
      if (npacked_per_int .eq. 4) then
         ind = 1
         do i = 1, nunpacked, 4
            ij = ior(ishft(unpacked(i  ),8),unpacked(i+1))
            kl = ior(ishft(unpacked(i+2),8),unpacked(i+3))
            ijkl = ior(ishft(ij,16),kl)
            packed(ind) = ijkl
            ind = ind + 1
         enddo
      else if (npacked_per_int .eq. 8) then
         ind = 1
         do i = 1, nunpacked, 8
            ij = ior(ishft(unpacked(i  ),8),unpacked(i+1))
            kl = ior(ishft(unpacked(i+2),8),unpacked(i+3))
            mn = ior(ishft(unpacked(i+4),8),unpacked(i+5))
            pq = ior(ishft(unpacked(i+6),8),unpacked(i+7))
            ijkl = ior(ishft(ij,16),kl)
            mnpq = ior(ishft(mn,16),pq)
            packed(ind) = ior(ishft(ijkl,32),mnpq)
            ind = ind + 1
         enddo
      else
         call errquit('util_pack_8: funny integer length',
     $        npacked_per_int, INPUT_ERR)
      endif
c
      end
      subroutine util_unpack_8(nunpacked, packed, unpacked)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
c
c     unpack 8 bit integers ... n must be a multiple of 8.
c
      integer nunpacked, packed(*), unpacked(*)
c
      integer npacked_per_int
      integer i, ijkl,ind,mnpq,mask
#include "bitops.fh"
c
      npacked_per_int = ma_sizeof(mt_int, 1, mt_byte)
c
c     The few extra iands take care of machines that insist on
c     sign extension when shifting
c
      if (npacked_per_int .eq. 4) then
         ind = 1
         do i = 1, nunpacked/4
            ijkl = packed(i)
            unpacked(ind  ) = iand(ishft(ijkl,-24),255)
            unpacked(ind+1) = iand(ishft(ijkl,-16),255)
            unpacked(ind+2) = iand(ishft(ijkl,-8),255)
            unpacked(ind+3) = iand(ijkl,255)
            ind = ind + 4
         enddo
      else if (npacked_per_int .eq. 8) then
         ind = 1
         mask = ishft(1,32)-1  ! 2^32 - 1
         do i = 1, nunpacked/8
            ijkl = iand(ishft(packed(i),-32),mask)
            mnpq = iand(packed(i),mask)
            unpacked(ind  ) = iand(ishft(ijkl,-24),255)
            unpacked(ind+1) = iand(ishft(ijkl,-16),255)
            unpacked(ind+2) = iand(ishft(ijkl, -8),255)
            unpacked(ind+3) = iand(ijkl,255)
            unpacked(ind+4) = iand(ishft(mnpq,-24),255)
            unpacked(ind+5) = iand(ishft(mnpq,-16),255)
            unpacked(ind+6) = iand(ishft(mnpq, -8),255)
            unpacked(ind+7) = iand(mnpq,255)
            ind = ind + 8
         enddo
      else
         call errquit('util_unpack_8 funny integer length',
     $        npacked_per_int, INPUT_ERR)
      endif
c
      end
      subroutine util_pack_16(nunpacked, packed, unpacked)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
c
c     pack 16 bit integers ... n must be a multiple of 4
c
      integer nunpacked, packed(*), unpacked(*)
      integer npacked_per_int
      integer i, ij, kl, ijkl, ind
#include "bitops.fh"
c
      npacked_per_int = ma_sizeof(mt_int, 1, mt_byte) / 2
c
      if (npacked_per_int .eq. 2) then
         ind = 1
         do i = 1, nunpacked, 2
            packed(ind) = ior(ishft(unpacked(i),16),unpacked(i+1))
            ind = ind + 1
         enddo
      else if (npacked_per_int .eq. 4) then
         ind = 1
         do i = 1, nunpacked, 4
            ij = ior(ishft(unpacked(i  ),16),unpacked(i+1))
            kl = ior(ishft(unpacked(i+2),16),unpacked(i+3))
            ijkl = ior(ishft(ij,32),kl)
            packed(ind) = ijkl
            ind = ind + 1
         enddo
      else
         call errquit('util_pack_16: funny integer length',
     $        npacked_per_int, INPUT_ERR)
      endif
c
      end
      subroutine util_unpack_16(nunpacked, packed, unpacked)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
c
c     unpack 16 bit integers ... n must be a multiple of 4.
c
      integer nunpacked, packed(*), unpacked(*)
c
      integer npacked_per_int
      integer i, ij, ijkl, ind
#include "bitops.fh"
c
      npacked_per_int = ma_sizeof(mt_int, 1, mt_byte) / 2
c
      if (npacked_per_int .eq. 2) then
         ind = 1
         do i = 1, nunpacked/2
            ij = packed(i)
            unpacked(ind)   = iand(ishft(ij,-16),65535)
            unpacked(ind+1) = iand(ij,65535)
            ind = ind + 2
         enddo
      else if (npacked_per_int .eq. 4) then
         ind = 1
         do i = 1, nunpacked/4
            ijkl = packed(i)
            unpacked(ind  ) = iand(ishft(ijkl,-48),65535)
            unpacked(ind+1) = iand(ishft(ijkl,-32),65535)
            unpacked(ind+2) = iand(ishft(ijkl,-16),65535)
            unpacked(ind+3) = iand(ijkl,65535)
            ind = ind + 4
         enddo
      else
         call errquit('util_unpack_16: funny integer length',
     $        npacked_per_int, INPUT_ERR)
      endif
c
      end
#endif
