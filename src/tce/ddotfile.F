      double precision function ddotfile(d_1,d_2,size)
c
c $Id$
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
cc      integer nxtask 
      integer next, nprocs, count
      INTEGER NXTASK
      EXTERNAL NXTASK
cc      external nxtask
      logical noloadbalance
c
c     new
c
      if ((fileisga(d_1).and.fileisga(d_2))) then
        ddotfile = ga_ddot(d_1,d_2)
        return
      endif
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
cc      next = nxtask(nprocs,1)
      next = NXTASK(nprocs, 1)
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
cc          next = nxtask(nprocs,1)
      next = NXTASK(nprocs, 1)
        endif
        count = count + 1
      enddo
cc      next = nxtask(-nprocs,1)
      next = NXTASK(-nprocs, 1)
      call ga_sync()
      if (.not.noloadbalance) call ga_dgop(mt_dbl,ddotfile,1,'+')
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
c
c
c
      double precision function ddotfile_1(d_1,d_2,size)
c
c $Id$
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
cc      integer nxtask 
      integer next, nprocs, count
      integer buffer_size_t1,size_t1
cc      external nxtask
      INTEGER NXTASK
      EXTERNAL NXTASK
      logical noloadbalance
c
c     new
c
      if ((fileisga(d_1).and.fileisga(d_2))) then
        ddotfile_1 = ga_ddot(d_1,d_2)
        return
      endif
c
      buffer_size_t1=tile_dim**2
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
      ddotfile_1 = 0.0d0
      nblocks = size/buffer_size_t1
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
cc      next = nxtask(nprocs,1)
      next = NXTASK(nprocs, 1)
      do i = 1,nblocks
        if (noloadbalance.or.(next.eq.count)) then
          if (util_print('ddotfile',print_debug))
     1      write(LuOut,9210) ga_nodeid(),(i-1)*blocksize,
     2                        min(size,i*blocksize)-(i-1)*blocksize
          call get_block(d_1,dbl_mb(k_temp1),
     1      min(size,i*blocksize)-(i-1)*blocksize,(i-1)*blocksize)
          call get_block(d_2,dbl_mb(k_temp2),
     1      min(size,i*blocksize)-(i-1)*blocksize,(i-1)*blocksize)
ccx          k = 0
ccx          do j = (i-1)*blocksize+1,min(size,i*blocksize)
ccx            k = k + 1
ccx            ddotfile_1 = ddotfile_1 + dbl_mb(k_temp1+k-1)
ccx     1                          * dbl_mb(k_temp2+k-1)
ccx          enddo
c
          size_t1=min(size,i*blocksize)-(i-1)*blocksize
          ddotfile_1=ddotfile_1
     1       +ddot(size_t1,dbl_mb(k_temp1),1,dbl_mb(k_temp2),1)
c
cc          next = nxtask(nprocs,1)
      next = NXTASK(nprocs, 1)
        endif
        count = count + 1
      enddo
cc      next = nxtask(-nprocs,1)
      next = NXTASK(-nprocs, 1)
      call ga_sync()
      if (.not.noloadbalance) call ga_dgop(mt_dbl,ddotfile_1,1,'+')
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
c
c
c
c
c
c
      double precision function ddotfile_2(d_1,d_2,size)
c
c $Id$
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
cc      integer nxtask
      integer  next, nprocs, count
      integer buffer_size_t2,size_t2
cc      external nxtask
      INTEGER NXTASK
      EXTERNAL NXTASK
      logical noloadbalance
c
c     new
c
      if ((fileisga(d_1).and.fileisga(d_2))) then
        ddotfile_2 = ga_ddot(d_1,d_2)
        return
      endif
c
      buffer_size_t2=tile_dim**4
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
      ddotfile_2 = 0.0d0
      nblocks = size/buffer_size_t2
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
cc      next = nxtask(nprocs,1)
      next = NXTASK(nprocs, 1)
      do i = 1,nblocks
        if (noloadbalance.or.(next.eq.count)) then
          if (util_print('ddotfile',print_debug))
     1      write(LuOut,9210) ga_nodeid(),(i-1)*blocksize,
     2                        min(size,i*blocksize)-(i-1)*blocksize
          call get_block(d_1,dbl_mb(k_temp1),
     1      min(size,i*blocksize)-(i-1)*blocksize,(i-1)*blocksize)
          call get_block(d_2,dbl_mb(k_temp2),
     1      min(size,i*blocksize)-(i-1)*blocksize,(i-1)*blocksize)
ccx          k = 0
ccx          do j = (i-1)*blocksize+1,min(size,i*blocksize)
ccx            k = k + 1
ccx            ddotfile_1 = ddotfile_1 + dbl_mb(k_temp1+k-1)
ccx     1                          * dbl_mb(k_temp2+k-1)
ccx          enddo
c
          size_t2=min(size,i*blocksize)-(i-1)*blocksize
          ddotfile_2=ddotfile_2
     1       +ddot(size_t2,dbl_mb(k_temp1),1,dbl_mb(k_temp2),1)
c
cc          next = nxtask(nprocs,1)
      next = NXTASK(nprocs, 1)
        endif
        count = count + 1
      enddo
cc      next = nxtask(-nprocs,1)
      next = NXTASK(-nprocs, 1)
      call ga_sync()
      if (.not.noloadbalance) call ga_dgop(mt_dbl,ddotfile_2,1,'+')
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
