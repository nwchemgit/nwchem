      subroutine add_block(d_file,array,size,offset)
c
c $Id: add_block.F,v 1.1 2002-10-18 23:42:05 sohirata Exp $
c
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "sf.fh"
#include "util.fh"
#include "stdio.fh"
#include "tce.fh"
      integer d_file
      integer size
      integer offset
      double precision array(size)
      integer l_temp,k_temp
      integer dummy
      integer request
      integer i
      logical nodezero
c
      nodezero = (ga_nodeid().eq.0)
      if (nodezero.and.util_print('add_block',print_debug)) then
         write(LuOut,9000) d_file,size,offset
      endif
      if (offset .eq. -1)
     1  call errquit('add_block: illegal offset',0)
      if (.not.ma_push_get(mt_dbl,size,'temporary',
     1  l_temp,k_temp))
     2  call errquit('add_block: MA problem',0)
      dummy=sf_read(d_file,dfloat(bytes*offset),
     1  dfloat(bytes*size),dbl_mb(k_temp),request)
      dummy=sf_wait(request)
      do i = 1,size
        dbl_mb(k_temp+i-1) = dbl_mb(k_temp+i-1) + array(i)
      enddo
      if (.not.ma_pop_stack(l_temp))
     1  call errquit('add_block: MA problem',0)
      dummy=sf_write(d_file,dfloat(bytes*offset),
     1  dfloat(bytes*size),dbl_mb(k_temp),request)
      dummy=sf_wait(request)
 9000 format(1x,'add_block request to file:',i10,
     1  ' size:',i10,' offset:',i10)
      return
      end
