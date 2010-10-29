      subroutine ma_pack(unpack,pack,n)
c
c $Id$
c
      implicit none
#include "mafdecls.fh"
#include "tce.fh"
#include "tce_main.fh"
      integer i,j,k,n
      integer g1b,g2b,ioff,joff
      double precision unpack(n,n)
      double precision pack(*)
      k = 0
      do g1b = 1,noa+nob+nva+nvb
        do g2b = g1b,noa+nob+nva+nvb
          joff = int_mb(k_offset+g1b-1)
          ioff = int_mb(k_offset+g2b-1)
          if (int_mb(k_spin+g1b-1).eq.int_mb(k_spin+g2b-1)) then
            do j = 1,int_mb(k_range+g1b-1)
              do i = 1,int_mb(k_range+g2b-1)
                k = k + 1
                pack(k) = unpack(j+joff,i+ioff)-unpack(i+ioff,j+joff)
              enddo
            enddo
          else
            do j = 1,int_mb(k_range+g1b-1)
              do i = 1,int_mb(k_range+g2b-1)
                k = k + 1
                pack(k) = unpack(j+joff,i+ioff)
              enddo
            enddo
          endif
        enddo
      enddo
      return
      end
