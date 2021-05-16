      subroutine daxpyfile(n,coefficients,d_in,d_out,size)
c
c $Id$
c
c     A linear combination of files added to the original
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
      integer n
      integer d_in(*)
      integer d_out
      integer size
      integer nblocks,blocksize
      integer l_temp1,k_temp1
      integer l_temp2,k_temp2
      integer i,j,k,l
      logical nodezero
      double precision coefficients(*)
      logical dinisga, doutisga
cc      integer nxtask 
      integer next, nprocs, count
cc      external nxtask
      INTEGER NXTASK
      EXTERNAL NXTASK
      logical noloadbalance
c
c     new
c
      if ( (fileisga(d_in(1)).and.fileisga(d_out)) ) then
        if (n.eq.1) then
c         d_c = alpha * d_a  +  beta * d_b
          call ga_add(coefficients(1),d_in(1), ! alpha,d_a
     1                1.0d0,d_out,             ! beta, d_b
     2                d_out)                   ! d_c
          return
        else
          do l = 1,n
            call ga_add(coefficients(l),d_in(l),1.0d0,d_out,d_out)
          enddo
          return
        endif
      endif
c
      if ((ioalg.eq.6).and.(fileisga(d_in(1)).neqv.fileisga(d_out)))
     1  call errquit('unable to daxpy mixed GA and EAF',0,UNKNOWN_ERR)
      noloadbalance = ((ioalg.eq.4).or.
     1                ((ioalg.eq.6).and.(.not.fileisga(d_in(1)))))
      dinisga = fileisga(d_in(1))
      doutisga = fileisga(d_out)
      nodezero = (ga_nodeid().eq.0)
      nblocks = size/buffer_size
      if (nblocks*buffer_size.lt.size) nblocks = nblocks + 1
      blocksize = size/nblocks
      if (nblocks*blocksize.lt.size) blocksize = blocksize + 1
      if (nodezero.and.util_print('daxpyfile',print_debug)) then
        do l = 1,n
          write(LuOut,9180) d_in(l),size,coefficients(l)
        enddo
        write(LuOut,9190) d_out,size
        write(LuOut,9200) nblocks,blocksize
      endif
      if (ioalg.eq.5) then
        if (.not.doutisga) call dratoga(d_out)
        do l = 1,n
          if (.not.dinisga) call dratoga(d_in(l))
          call ga_add(coefficients(l),filehandles(d_in(l)),1.0d0,
     1         filehandles(d_out),filehandles(d_out))
          if (.not.dinisga) call gatodra(d_in(l))
        enddo
        if (.not.doutisga) call gatodra(d_out)
      else
        if (.not.ma_push_get(mt_dbl,blocksize,'temporary',
     1    l_temp1,k_temp1))
     2    call errquit('daxpyfile: MA problem',0,MA_ERR)
        if (.not.ma_push_get(mt_dbl,blocksize,'temporary',
     1    l_temp2,k_temp2))
     2    call errquit('daxpyfile: MA problem',1,MA_ERR)
        nprocs = GA_NNODES()
        count = 0
cc        next = nxtask(nprocs,1)
      next = NXTASK(nprocs, 1)
        do i = 1,nblocks
          if (noloadbalance.or.(next.eq.count)) then
            if (util_print('daxpyfile',print_debug))
     1        write(LuOut,9210) ga_nodeid(),(i-1)*blocksize,
     2                          min(size,i*blocksize)-(i-1)*blocksize
            k = 0
            do j = (i-1)*blocksize+1,min(size,i*blocksize)
              k = k + 1
              dbl_mb(k_temp1+k-1) = 0.0d0
            enddo
            do l = 1,n
              call get_block(d_in(l),dbl_mb(k_temp2),
     1          min(size,i*blocksize)-(i-1)*blocksize,(i-1)*blocksize)
              k = 0
              do j = (i-1)*blocksize+1,min(size,i*blocksize)
                k = k + 1
                dbl_mb(k_temp1+k-1) = dbl_mb(k_temp1+k-1) 
     1          + coefficients(l) * dbl_mb(k_temp2+k-1)
              enddo
            enddo
            call add_block(d_out,dbl_mb(k_temp1),
     1        min(size,i*blocksize)-(i-1)*blocksize,(i-1)*blocksize)
cc            next = nxtask(nprocs,1)
      next = NXTASK(nprocs, 1)
          endif
          count = count + 1
        enddo
cc        next = nxtask(-nprocs,1)
      next = NXTASK(-nprocs, 1)
        if (.not.ma_pop_stack(l_temp2))
     1    call errquit('createfile: MA problem',2,MA_ERR)
        if (.not.ma_pop_stack(l_temp1))
     1    call errquit('createfile: MA problem',3,MA_ERR)
      endif
      call ga_sync()
 9180 format(/,1x,'  input file handle:',i10,' size:',i10,
     1            ' coefficient:',f15.10)
 9190 format(/,1x,' output file handle:',i10,' size:',i10)
 9200 format(  1x,'            nblocks:',i10,' size:',i10)
 9210 format(  1x,i3,'          offset:',i10,' size:',i10)
      return
      end
