      subroutine util_pack_16(nunpacked, packed, unpacked)
      implicit none
#include "mafdecls.fh"
c
c     pack 16 bit integers ... n must be a multiple of 4
c
      integer nunpacked, packed(*), unpacked(*)
      integer npacked_per_int
      integer i, ij, kl, ijkl, ind
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
     $        npacked_per_int)
      endif
c
      end
      subroutine util_unpack_16(nunpacked, packed, unpacked)
      implicit none
#include "mafdecls.fh"
c
c     unpack 16 bit integers ... n must be a multiple of 4.
c
      integer nunpacked, packed(*), unpacked(*)
c
      integer npacked_per_int
      integer i, ij, ijkl, ind
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
     $        npacked_per_int)
      endif
c
      end
      subroutine util_pack_8(nunpacked, packed, unpacked)
      implicit none
#include "mafdecls.fh"
c
c     pack 8 bit integers ... n must be a multiple of 8
c
      integer nunpacked, packed(*), unpacked(*)
      integer npacked_per_int
      integer i, ij, kl, ijkl, ind, mn, pq, mnpq
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
     $        npacked_per_int)
      endif
c
      end
      subroutine util_unpack_8(nunpacked, packed, unpacked)
      implicit none
#include "mafdecls.fh"
c
c     unpack 8 bit integers ... n must be a multiple of 8.
c
      integer nunpacked, packed(*), unpacked(*)
c
      integer npacked_per_int
      integer i, ijkl,ind,mnpq,mask
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
*            unpacked(ind  ) = iand(ishft(ijkl,-24),255)
            unpacked(ind  ) =      ishft(ijkl,-24)
            unpacked(ind+1) = iand(ishft(ijkl,-16),255)
            unpacked(ind+2) = iand(ishft(ijkl, -8),255)
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
     $        npacked_per_int)
      endif
c
      end
