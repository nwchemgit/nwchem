      subroutine tce_zero(d_sf,size)
c
c $Id$
c
c     Zero scratch a file
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
      integer d_sf
      integer size
      integer nblocks,blocksize
      integer l_temp,k_temp
      integer i,j,k
      logical nodezero
cc      integer nxtask 
      integer next, nprocs, count
cc      external nxtask
      INTEGER NXTASK
      EXTERNAL NXTASK
      logical noloadbalance
c
      nodezero = (ga_nodeid().eq.0)
      noloadbalance = ((ioalg.eq.4).or.
     1                ((ioalg.eq.6).and.(.not.fileisga(d_sf))))
      nblocks = size/buffer_size
      if (nblocks*buffer_size.lt.size) nblocks = nblocks + 1
      blocksize = size/nblocks
      if (nblocks*blocksize.lt.size) blocksize = blocksize + 1
      if (nodezero.and.util_print('zero',print_debug)) then
        write(LuOut,9190) d_sf,size
        write(LuOut,9200) nblocks,blocksize
      endif
      if (ioalg.eq.2) then
        call ga_zero(d_sf)
      else if (ioalg.eq.5) then
        if (fileisga(d_sf)) then
          call ga_zero(filehandles(d_sf))
        else
          call dratoga(d_sf)
          call ga_zero(filehandles(d_sf))
          call gatodra(d_sf)
        endif
      else if ((ioalg.eq.6).and.fileisga(d_sf)) then
        call ga_zero(filehandles(d_sf))
      else
        if (.not.ma_push_get(mt_dbl,blocksize,'temporary',
     1    l_temp,k_temp))
     2    call errquit('tce_zero: MA problem',0,MA_ERR)
        nprocs = GA_NNODES()
        count = 0
cc        next = nxtask(nprocs,1)
      next = NXTASK(nprocs, 1)
        do i = 1,nblocks
          if (noloadbalance.or.(next.eq.count)) then
            k = 0
            do j = (i-1)*blocksize+1,min(size,i*blocksize)
              k = k + 1
              dbl_mb(k_temp+k-1) = 0.0d0
            enddo
            if (util_print('zero',print_debug))
     1        write(LuOut,9210) ga_nodeid(),(i-1)*blocksize,
     2                          min(size,i*blocksize)-(i-1)*blocksize
            call put_block(d_sf,dbl_mb(k_temp),
     1        min(size,i*blocksize)-(i-1)*blocksize,(i-1)*blocksize)
cc            next = nxtask(nprocs,1)
      next = NXTASK(nprocs, 1)
          endif
          count = count + 1
        enddo
cc        next = nxtask(-nprocs,1)
      next = NXTASK(-nprocs, 1)
        if (.not.ma_pop_stack(l_temp))
     1    call errquit('createfile: MA problem',1,MA_ERR)
        call ga_sync()
      endif
 9190 format(/,1x,'Zero scratch handle:',i10,' size:',i10)
 9200 format(  1x,'            nblocks:',i10,' size:',i10)
 9210 format(  1x,i3,'          offset:',i10,' size:',i10)
      return
      end















      subroutine tce_zero_sf(d_sf,size)
c
c $Id$
c
c     Zero scratch a file
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
      integer d_sf
      integer size
      integer nblocks,blocksize
      integer l_temp,k_temp
      integer i,j,k
      logical nodezero
cc      integer nxtask 
      integer next, nprocs, count
cc      external nxtask
      INTEGER NXTASK
      EXTERNAL NXTASK
      logical noloadbalance
c
      nodezero = (ga_nodeid().eq.0)
      noloadbalance = .false.
      nblocks = size/buffer_size
      if (nblocks*buffer_size.lt.size) nblocks = nblocks + 1
      blocksize = size/nblocks
      if (nblocks*blocksize.lt.size) blocksize = blocksize + 1
      if (nodezero.and.util_print('zero',print_debug)) then
        write(LuOut,9190) d_sf,size
        write(LuOut,9200) nblocks,blocksize
      endif
c
c      else
c
        if (.not.ma_push_get(mt_dbl,blocksize,'temporary',
     1    l_temp,k_temp))
     2    call errquit('tce_zero: MA problem',0,MA_ERR)
        nprocs = GA_NNODES()
        count = 0
cc        next = nxtask(nprocs,1)
        next = NXTASK(nprocs, 1)
        do i = 1,nblocks
          if (noloadbalance.or.(next.eq.count)) then
            k = 0
            do j = (i-1)*blocksize+1,min(size,i*blocksize)
              k = k + 1
              dbl_mb(k_temp+k-1) = 0.0d0
            enddo
            if (util_print('zero',print_debug))
     1        write(LuOut,9210) ga_nodeid(),(i-1)*blocksize,
     2                          min(size,i*blocksize)-(i-1)*blocksize
c
c            call put_block(d_sf,dbl_mb(k_temp),
c
            call put_block_sf(d_sf,dbl_mb(k_temp),
c
     1        min(size,i*blocksize)-(i-1)*blocksize,(i-1)*blocksize)
cc            next = nxtask(nprocs,1)
            next = NXTASK(nprocs, 1)
          endif
          count = count + 1
        enddo
cc        next = nxtask(-nprocs,1)
        next = NXTASK(-nprocs, 1)
        if (.not.ma_pop_stack(l_temp))
     1    call errquit('createfile: MA problem',1,MA_ERR)
        call ga_sync()
c
c      endif
c
 9190 format(/,1x,'Zero scratch handle:',i10,' size:',i10)
 9200 format(  1x,'            nblocks:',i10,' size:',i10)
 9210 format(  1x,i3,'          offset:',i10,' size:',i10)
      return
      end
