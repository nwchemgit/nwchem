      subroutine sf_print(sf_handle,size)
c
c $Id: sf_print.F,v 1.3 2002-10-29 00:58:06 sohirata Exp $
c
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "util.fh"
#include "sf.fh"
#include "stdio.fh"
#include "tce.fh"
      integer sf_handle
      integer size
      integer request
      integer l_sf,k_sf
      integer i
      double precision checksum
c
      call ga_sync()
      call ga_lock(0)
      if (.not.ma_push_get(mt_dbl,size,'sf',l_sf,k_sf))
     1  call errquit('sf_print: MA problem',0)
      if (sf_read(sf_handle,0.0d0,dfloat(bytes*size),
     1  dbl_mb(k_sf),request).ne.0) 
     2  call errquit('sf_print: sf problem',0)
      if (sf_wait(request).ne.0)
     1  call errquit('sf_print: sf problem',0)
      call ga_unlock(0)
      if (ga_nodeid().eq.0) then
        call ma_print(dbl_mb(k_sf),size,1,'Shared file')
        checksum = 0.0d0
        do i = 1,size
          checksum = checksum + dbl_mb(k_sf+i-1)**2
        enddo
        write(LuOut,9000) checksum
      endif
      if (.not.ma_pop_stack(l_sf))
     1  call errquit('sf_print: MA problem',0)
 9000 format(1x,'Checksum = ',f15.10)
      return
      end
