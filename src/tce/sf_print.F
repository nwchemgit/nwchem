      subroutine sf_print(sf_handle,size)
c
c $Id: sf_print.F,v 1.1 2002-10-18 23:42:07 sohirata Exp $
c
      implicit none
#include "mafdecls.fh"
#include "util.fh"
#include "sf.fh"
#include "tce.fh"
      integer sf_handle
      integer size
      integer request
      integer dummy
      integer l_sf,k_sf
c
      if (.not.ma_push_get(mt_dbl,size,'sf',l_sf,k_sf))
     1  call errquit('sf_print: MA problem',0)
      dummy = sf_read(sf_handle,0.0d0,dfloat(bytes*size),
     1  dbl_mb(k_sf),request)
      dummy = sf_wait(request)
      call ma_print(dbl_mb(k_sf),size,1,'Shared file')
      if (.not.ma_pop_stack(l_sf))
     1  call errquit('sf_print: MA problem',0)
      return
      end
