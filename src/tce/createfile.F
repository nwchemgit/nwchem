      subroutine createfile(filename,d_sf,size)
c
c $Id: createfile.F,v 1.2 2002-10-23 01:38:51 sohirata Exp $
c
c     Creates a shared file and returns the handle
c
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "util.fh"
#include "sf.fh"
#include "tce.fh"
      character*(*) filename
      integer d_sf
      integer size
      integer l_temp,k_temp
      integer i
c
      if (sf_create(filename,dfloat(size),dfloat(size),
     1  1.0d4,d_sf).ne.0)
     2  call errquit('createfile: sf problem',0)
      if (.not.ma_push_get(mt_dbl,size,'temporary',
     1  l_temp,k_temp))
     2  call errquit('createfile: MA problem',0)
      do i = 1,size
        dbl_mb(k_temp+i-1) = 0.0d0
      enddo
      call put_block(d_sf,dbl_mb(k_temp),size,0)
      if (.not.ma_pop_stack(l_temp))
     1  call errquit('createfile: MA problem',1)
      call ga_sync()
      return
      end
