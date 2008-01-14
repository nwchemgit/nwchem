      double precision function ddotfile(d_1,d_2,size)
c
c $Id: ddotfile.F,v 1.4 2008-01-14 15:20:32 jhammond Exp $
c
c     Vector inner product of two files of the same size
c
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "util.fh"
#include "sf.fh"
#include "msgids.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
      integer d_1
      integer d_2
      integer size
      integer nblocks,blocksize
      integer l_temp1,k_temp1
      integer l_temp2,k_temp2
      integer i,j,k
      logical nodezero
      logical d_1isga, d_2isga
      integer nxtval, next, nprocs, count
      external nxtval
      logical noloadbalance
c
      if ((ioalg.eq.6).and.(fileisga(d_1).neqv.fileisga(d_2)))
     1  call errquit('unable to ddot GA and EAF',0,UNKNOWN_ERR)
      noloadbalance = ((ioalg.eq.4).or.
     1                ((ioalg.eq.6).and.(.not.fileisga(d_1))))
      d_1isga = fileisga(d_1)
      d_2isga = fileisga(d_2)
      if (.not.d_1isga) call dratoga(d_1)
      if ((.not.d_2isga).and.(d_1.ne.d_2)) call dratoga(d_2)
c     if ((ioalg.eq.2).or.((ioalg.eq.5).and.fileisga(d_1))
c    1                .or.((ioalg.eq.6).and.fileisga(d_1))) then
c       ddotfile = ga_ddot(d_1,d_2)
c       if (.not.d_1isga) call gatodra(d_1)
c       if ((.not.d_2isga).and.(d_1.ne.d_2)) call gatodra(d_2)
c       return
c     endif
      nodezero = (ga_nodeid().eq.0)
      ddotfile = 0.0d0
      nblocks = size/buffer_size
      if (nblocks*buffer_size.lt.size) nblocks = nblocks + 1
      blocksize = size/nblocks
      if (nblocks*blocksize.lt.size) blocksize = blocksize + 1
      if (nodezero.and.util_print('ddotfile',print_debug)) then
        write(LuOut,9190) d_1,size
        write(LuOut,9190) d_2,size
        write(LuOut,9200) nblocks,blocksize
      endif
      if (.not.ma_push_get(mt_dbl,blocksize,'temporary',
     1  l_temp1,k_temp1))
     2  call errquit('ddotfile: MA problem',0,MA_ERR)
      if (.not.ma_push_get(mt_dbl,blocksize,'temporary',
     1  l_temp2,k_temp2))
     2  call errquit('ddotfile: MA problem',1,MA_ERR)
      nprocs = GA_NNODES()
      count = 0
      next = NXTVAL(nprocs)
      do i = 1,nblocks
        if (noloadbalance.or.(next.eq.count)) then
          if (util_print('ddotfile',print_debug))
     1      write(LuOut,9210) ga_nodeid(),(i-1)*blocksize,
     2                        min(size,i*blocksize)-(i-1)*blocksize
          call get_block(d_1,dbl_mb(k_temp1),
     1      min(size,i*blocksize)-(i-1)*blocksize,(i-1)*blocksize)
          call get_block(d_2,dbl_mb(k_temp2),
     1      min(size,i*blocksize)-(i-1)*blocksize,(i-1)*blocksize)
          k = 0
          do j = (i-1)*blocksize+1,min(size,i*blocksize)
            k = k + 1
            ddotfile = ddotfile + dbl_mb(k_temp1+k-1)
     1                          * dbl_mb(k_temp2+k-1)
          enddo
          next = NXTVAL(nprocs)
        endif
        count = count + 1
      enddo
      next = NXTVAL(-nprocs)
      call ga_sync()
      if (.not.noloadbalance) call ga_dgop(msg_tcemem,ddotfile,1,'+')
      if (.not.d_1isga) call gatodra(d_1)
      if ((.not.d_2isga).and.(d_1.ne.d_2)) call gatodra(d_2)
      if (.not.ma_pop_stack(l_temp2))
     1  call errquit('createfile: MA problem',2,MA_ERR)
      if (.not.ma_pop_stack(l_temp1))
     1  call errquit('createfile: MA problem',3,MA_ERR)
 9190 format(/,1x,'        file handle:',i10,' size:',i10)
 9200 format(  1x,'            nblocks:',i10,' size:',i10)
 9210 format(  1x,i3,'          offset:',i10,' size:',i10)
      return
      end
