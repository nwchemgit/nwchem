      subroutine createfile(filename,d_sf,size)
c
c $Id: createfile.F,v 1.1 2002-10-18 23:42:06 sohirata Exp $
c
c     Creates a shared file and returns the handle
c
      implicit none
#include "mafdecls.fh"
#include "stdio.fh"
#include "util.fh"
#include "sf.fh"
#include "tce.fh"
      character*(*) filename
      integer d_sf
      integer size
      integer dummy
      integer l_temp,k_temp
      integer i
c
      dummy = sf_create(filename,dfloat(size),dfloat(size),
     1  1.0d4,d_sf)
      if (.not.ma_push_get(mt_dbl,size,'temporary',
     1  l_temp,k_temp))
     2  call errquit('createfile: MA problem',0)
      do i = 1,size
        dbl_mb(k_temp+i-1) = 0.0d0
      enddo
      call put_block(d_sf,dbl_mb(k_temp),size,0)
      if (.not.ma_pop_stack(l_temp))
     1  call errquit('createfile: MA problem',1)
      return
      end
