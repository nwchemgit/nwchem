      subroutine copyfile(d_1,d_2,size)
c
c $Id$
c
c     Duplicate a file
c
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "util.fh"
#include "sf.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
      integer d_1
      integer d_2
      integer size
      integer nblocks,blocksize
      integer l_temp,k_temp
      integer i
      logical nodezero
      logical d_1isga
      logical d_2isga
      integer  next, nprocs, count
      INTEGER NXTASK
      EXTERNAL NXTASK
      logical noloadbalance
c
c     new
c
      if ( fileisga(d_1).and.fileisga(d_2) ) then
        call ga_copy(d_1,d_2)
        return
      endif
c
      if ((ioalg.eq.6).and.(fileisga(d_1).neqv.fileisga(d_2)))
     1  call errquit('unable to copy GA to/from EAF',0,UNKNOWN_ERR)
      noloadbalance = ((ioalg.eq.4).or.
     1                ((ioalg.eq.6).and.(.not.fileisga(d_1))))
      d_1isga = fileisga(d_1)
      d_2isga = fileisga(d_2)
      if (.not.d_1isga) call dratoga(d_1)
      if (.not.d_2isga) call dratoga(d_2)
      nodezero = (ga_nodeid().eq.0)
      nblocks = size/buffer_size
      if (nblocks*buffer_size.lt.size) nblocks = nblocks + 1
      blocksize = size/nblocks
      if (nblocks*blocksize.lt.size) blocksize = blocksize + 1
      if (nodezero.and.util_print('copyfile',print_debug)) then
        write(LuOut,9190) d_1,size
        write(LuOut,9190) d_2,size
        write(LuOut,9200) nblocks,blocksize
      endif
      if (.not.ma_push_get(mt_dbl,blocksize,'temporary',
     1  l_temp,k_temp))
     2  call errquit('copyfile: MA problem',0,MA_ERR)
      nprocs = GA_NNODES()
      count = 0
cc      next = nxtask(nprocs,1)
      next = NXTASK(nprocs, 1)
      do i = 1,nblocks
        if (noloadbalance.or.(next.eq.count)) then
          if (util_print('copyfile',print_debug))
     1      write(LuOut,9210) ga_nodeid(),(i-1)*blocksize,
     2                        min(size,i*blocksize)-(i-1)*blocksize
          call get_block(d_1,dbl_mb(k_temp),
     1      min(size,i*blocksize)-(i-1)*blocksize,(i-1)*blocksize)
          call put_block(d_2,dbl_mb(k_temp),
     1      min(size,i*blocksize)-(i-1)*blocksize,(i-1)*blocksize)
cc          next = nxtask(nprocs,1)
      next = NXTASK(nprocs, 1)
        endif
        count = count + 1
      enddo
cc      next = nxtask(-nprocs,1)
      next = NXTASK(-nprocs, 1)
      if (.not.ma_pop_stack(l_temp))
     1  call errquit('copyfile: MA problem',1,MA_ERR)
      call ga_sync()
      if (.not.d_1isga) call gatodra(d_1)
      if (.not.d_2isga) call gatodra(d_2)
 9190 format(/,1x,'        file handle:',i10,' size:',i10)
 9200 format(  1x,'            nblocks:',i10,' size:',i10)
 9210 format(  1x,i3,'          offset:',i10,' size:',i10)
      return
      end
