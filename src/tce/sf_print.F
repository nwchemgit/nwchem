      subroutine sf_print(sf_handle,size)
c
c $Id: sf_print.F,v 1.2 2002-10-23 01:38:52 sohirata Exp $
c
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "util.fh"
#include "sf.fh"
#include "tce.fh"
      integer sf_handle
      integer size
      integer request
      integer l_sf,k_sf
c
      call ga_sync()
      if (ga_nodeid().eq.0) then
        if (.not.ma_push_get(mt_dbl,size,'sf',l_sf,k_sf))
     1    call errquit('sf_print: MA problem',0)
        if (sf_read(sf_handle,0.0d0,dfloat(bytes*size),
     1    dbl_mb(k_sf),request).ne.0) 
     2    call errquit('sf_print: sf problem',0)
        if (sf_wait(request).ne.0)
     1    call errquit('sf_print: sf problem',0)
        call ma_print(dbl_mb(k_sf),size,1,'Shared file')
        if (.not.ma_pop_stack(l_sf))
     1    call errquit('sf_print: MA problem',0)
      endif
      return
      end
