      subroutine copyfile(d_1,d_2,size)
c
c $Id: copyfile.F,v 1.5 2003-09-16 21:50:11 sohirata Exp $
c
c     Duplicate a file
c
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "util.fh"
#include "sf.fh"
#include "tce.fh"
      integer d_1
      integer d_2
      integer size
      integer nblocks,blocksize
      integer l_temp,k_temp
      integer i
      logical nodezero
c
      nodezero = (ga_nodeid().eq.0)
      nblocks = size/buffer_size + 1
      blocksize = size/nblocks + 1
      if (nodezero.and.util_print('copyfile',print_debug)) then
        write(LuOut,9190) d_1,size
        write(LuOut,9190) d_2,size
        write(LuOut,9200) nblocks,blocksize
      endif
      if (.not.ma_push_get(mt_dbl,blocksize,'temporary',
     1  l_temp,k_temp))
     2  call errquit('copyfile: MA problem',0)
      do i = 1,nblocks
        if (nodezero.and.util_print('copyfile',print_debug))
     1    write(LuOut,9210) (i-1)*blocksize,
     2                      min(size,i*blocksize)-(i-1)*blocksize
        call get_block(d_1,dbl_mb(k_temp),
     1    min(size,i*blocksize)-(i-1)*blocksize,(i-1)*blocksize)
        call put_block(d_2,dbl_mb(k_temp),
     1    min(size,i*blocksize)-(i-1)*blocksize,(i-1)*blocksize)
      enddo
      if (.not.ma_pop_stack(l_temp))
     1  call errquit('copyfile: MA problem',1)
      call ga_sync()
 9190 format(/,1x,'        file handle:',i10,' size:',i10)
 9200 format(1x,'            nblocks:',i10,' size:',i10)
 9210 format(1x,'             offset:',i10,' size:',i10)
      return
      end
