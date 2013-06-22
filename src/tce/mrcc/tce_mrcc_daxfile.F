      subroutine mrcc_daxfile(n,coefficients,d_in,d_out,size,iref)
c
c     A linear combination of files
c
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "util.fh"
#include "sf.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_mrcc.fh"
#include "tce_mrcc_diis.fh"
#include "tce_main.fh"

      integer n
      integer d_in(maxdiism,maxref)
      integer d_out
      integer size
      integer nblocks,blocksize
      integer l_temp1,k_temp1
      integer l_temp2,k_temp2
      integer i,j,k,l
      logical nodezero
      double precision coefficients(*)
      logical dinisga, doutisga
      integer nxtask, next, nprocs, count
      external nxtask
      integer nxtasksub
      external nxtasksub
      logical noloadbalance
      integer iref

      if (ioalg.ne.2)
     1  call errquit('unable to dax non-GA file',0,UNKNOWN_ERR)

      nodezero = (ga_nodeid().eq.0)

      nblocks = size/buffer_size

      if (nblocks*buffer_size.lt.size) nblocks = nblocks + 1
      blocksize = size/nblocks
      if (nblocks*blocksize.lt.size) blocksize = blocksize + 1

cjb HERE
        if (.not.ma_push_get(mt_dbl,blocksize,'temporary',
     1    l_temp1,k_temp1))
     2    call errquit('daxfile: MA problem',0,MA_ERR)

        if (.not.ma_push_get(mt_dbl,blocksize,'temporary',
     1    l_temp2,k_temp2))
     2    call errquit('daxfile: MA problem',1,MA_ERR)

        count = 0
        if(lusesub) then
          nprocs = GA_pgroup_NNODES(mypgid)
          next = NXTASKsub(nprocs,1,mypgid)
        else
          nprocs = GA_NNODES()
          next = NXTASK(nprocs,1)
        endif

        do i = 1,nblocks
          if (next.eq.count) then
            if (util_print('daxfile',print_debug))
     1        write(LuOut,9210) ga_nodeid(),(i-1)*blocksize,
     2                          min(size,i*blocksize)-(i-1)*blocksize
            k = 0
            do j = (i-1)*blocksize+1,min(size,i*blocksize)
              k = k + 1
              dbl_mb(k_temp1+k-1) = 0.0d0
            enddo
            do l = 1,n
              call get_block(d_in(l,iref),dbl_mb(k_temp2),
     1          min(size,i*blocksize)-(i-1)*blocksize,(i-1)*blocksize)
              k = 0
              do j = (i-1)*blocksize+1,min(size,i*blocksize)
                k = k + 1
                dbl_mb(k_temp1+k-1) = dbl_mb(k_temp1+k-1) 
     1          + coefficients(l) * dbl_mb(k_temp2+k-1)
              enddo
            enddo
            call put_block(d_out,dbl_mb(k_temp1),
     1        min(size,i*blocksize)-(i-1)*blocksize,(i-1)*blocksize)
          if(lusesub) then
            next = NXTASKsub(nprocs,1,mypgid)
          else
            next = NXTASK(nprocs,1)
          endif
          endif
          count = count + 1
        enddo

          if(lusesub) then
            next = NXTASKsub(-nprocs,1,mypgid)
            call ga_pgroup_sync(mypgid)
          else
            next = NXTASK(-nprocs,1)
            call ga_sync()
          endif

        if (.not.ma_pop_stack(l_temp2))
     1    call errquit('createfile: MA problem',2,MA_ERR)

        if (.not.ma_pop_stack(l_temp1))
     1    call errquit('createfile: MA problem',3,MA_ERR)

 9180 format(/,1x,'  input file handle:',i10,' size:',i10,
     1            ' coefficient:',f15.10)
 9190 format(/,1x,' output file handle:',i10,' size:',i10)
 9200 format(  1x,'            nblocks:',i10,' size:',i10)
 9210 format(  1x,i3,'          offset:',i10,' size:',i10)
      return
      end
c $Id$
