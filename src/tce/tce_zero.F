      subroutine tce_zero(d_sf,size)
c
c $Id: tce_zero.F,v 1.5 2002-12-17 00:03:10 sohirata Exp $
c
c     Zero scratch a file
c
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "util.fh"
#include "sf.fh"
#include "tce.fh"
      integer d_sf
      integer size
      integer nblocks,blocksize
      integer l_temp,k_temp
      integer i,j,k
      logical nodezero
c
      nodezero = (ga_nodeid().eq.0)
      nblocks = size/buffer_size + 1
      blocksize = size/nblocks + 1
      if (nodezero.and.util_print('zero',print_debug)) then
        write(LuOut,9190) d_sf,size
        write(LuOut,9200) nblocks,blocksize
      endif
      if (.not.ma_push_get(mt_dbl,blocksize,'temporary',
     1  l_temp,k_temp))
     2  call errquit('tce_zero: MA problem',0)
      do i = 1,nblocks
        k = 0
        do j = (i-1)*blocksize+1,min(size,i*blocksize)
          k = k + 1
          dbl_mb(k_temp+k-1) = 0.0d0
        enddo
        if (nodezero.and.util_print('zero',print_debug))
     1    write(LuOut,9210) (i-1)*blocksize,
     2                      min(size,i*blocksize)-(i-1)*blocksize
        call put_block(d_sf,dbl_mb(k_temp),
     1    min(size,i*blocksize)-(i-1)*blocksize,(i-1)*blocksize)
      enddo
      if (.not.ma_pop_stack(l_temp))
     1  call errquit('createfile: MA problem',1)
      call ga_sync()
 9190 format(/,1x,'Zero scratch handle:',i10,' size:',i10)
 9200 format(1x,'            nblocks:',i10,' size:',i10)
 9210 format(1x,'             offset:',i10,' size:',i10)
      return
      end
