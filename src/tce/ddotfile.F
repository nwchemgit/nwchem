      double precision function ddotfile(d_1,d_2,size)
c
c $Id: ddotfile.F,v 1.1 2002-12-24 06:13:31 sohirata Exp $
c
c     Vector inner product of two files of the same size
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
      integer l_temp1,k_temp1
      integer l_temp2,k_temp2
      integer i,j,k
      logical nodezero
c
      nodezero = (ga_nodeid().eq.0)
      ddotfile = 0.0d0
      nblocks = size/buffer_size + 1
      blocksize = size/nblocks + 1
      if (nodezero.and.util_print('ddotfile',print_debug)) then
        write(LuOut,9190) d_1,size
        write(LuOut,9190) d_2,size
        write(LuOut,9200) nblocks,blocksize
      endif
      if (.not.ma_push_get(mt_dbl,blocksize,'temporary',
     1  l_temp1,k_temp1))
     2  call errquit('tce_zero: MA problem',0)
      if (.not.ma_push_get(mt_dbl,blocksize,'temporary',
     1  l_temp2,k_temp2))
     2  call errquit('tce_zero: MA problem',1)
      do i = 1,nblocks
        if (nodezero.and.util_print('ddotfile',print_debug))
     1    write(LuOut,9210) (i-1)*blocksize,
     2                      min(size,i*blocksize)-(i-1)*blocksize
        call get_block(d_1,dbl_mb(k_temp1),
     1    min(size,i*blocksize)-(i-1)*blocksize,(i-1)*blocksize)
        call get_block(d_2,dbl_mb(k_temp2),
     1    min(size,i*blocksize)-(i-1)*blocksize,(i-1)*blocksize)
        k = 0
        do j = (i-1)*blocksize+1,min(size,i*blocksize)
          k = k + 1
          ddotfile = ddotfile + dbl_mb(k_temp1+k-1)*dbl_mb(k_temp2+k-1)
        enddo
      enddo
      if (.not.ma_pop_stack(l_temp2))
     1  call errquit('createfile: MA problem',2)
      if (.not.ma_pop_stack(l_temp1))
     1  call errquit('createfile: MA problem',3)
      call ga_sync()
 9190 format(/,1x,'        file handle:',i10,' size:',i10)
 9200 format(1x,'            nblocks:',i10,' size:',i10)
 9210 format(1x,'             offset:',i10,' size:',i10)
      return
      end
