      subroutine sf_print(sf_handle,size)
c
c $Id: sf_print.F,v 1.4 2002-10-30 21:49:14 sohirata Exp $
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
      if (.not.ma_push_get(mt_dbl,size,'sf',l_sf,k_sf))
     1  call errquit('sf_print: MA problem',0)
      call get_block(sf_handle,dbl_mb(k_sf),size,0)
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
