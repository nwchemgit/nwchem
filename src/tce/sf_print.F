      subroutine sf_print(sf_handle,size)
c
c $Id$
c
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "util.fh"
#include "sf.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
      integer sf_handle
      integer size
      integer l_sf,k_sf
      integer i
      logical wasdra
      double precision checksum
c
      if (.not.ma_push_get(mt_dbl,size,'sf',l_sf,k_sf))
     1  call errquit('sf_print: MA problem',0,DISK_ERR)
      if (.not.fileisga(sf_handle)) then
        call dratoga(sf_handle)
        wasdra = .true.
      else
        wasdra = .false.
      endif
      call get_block(sf_handle,dbl_mb(k_sf),size,0)
      if (wasdra) call gatodra(sf_handle)
      write(LuOut,9000) ga_nodeid()
c     call ma_print(dbl_mb(k_sf),size,1,'Global file')
      write(LuOut,9020)
      do i = 1,size
        if (dabs(dbl_mb(k_sf+i-1)).gt.1.0d-10)
     1  write(LuOut,9030) dbl_mb(k_sf+i-1)
      enddo
      write(LuOut,9020)
      checksum = 0.0d0
      do i = 1,size
        checksum = checksum + dbl_mb(k_sf+i-1)**2
      enddo
      write(LuOut,9010) checksum
      if (.not.ma_pop_stack(l_sf))
     1  call errquit('sf_print: MA problem',0,MA_ERR)
 9000 format(/,1x,'Node = ',i3)
 9010 format(1x,'Checksum = ',f15.10)
 9020 format(1x,'---------------')
 9030 format(1x,f15.10)
      return
      end
