      subroutine daxpyfile(n,coefficients,d_in,d_out,size)
c
c $Id: daxpyfile.F,v 1.4 2003-08-05 00:37:51 sohirata Exp $
c
c     A linear combination of files added to the original
c
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "util.fh"
#include "sf.fh"
#include "tce.fh"
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
c
      nodezero = (ga_nodeid().eq.0)
      nblocks = size/buffer_size + 1
      blocksize = size/nblocks + 1
      if (nodezero.and.util_print('daxpyfile',print_debug)) then
        do l = 1,n
          write(LuOut,9180) d_in(l),size,coefficients(l)
        enddo
        write(LuOut,9190) d_out,size
        write(LuOut,9200) nblocks,blocksize
      endif
      if (ioalg.eq.5) then
        call dratoga(d_out)
        do l = 1,n
          call dratoga(d_in(l))
          call ga_add(coefficients(l),filehandles(d_in(l)),1.0d0,
     1         filehandles(d_out),filehandles(d_out))
          call gatodra(d_in(l))
        enddo
        call gatodra(d_out)
      else
        if (.not.ma_push_get(mt_dbl,blocksize,'temporary',
     1    l_temp1,k_temp1))
     2    call errquit('daxpyfile: MA problem',0)
        if (.not.ma_push_get(mt_dbl,blocksize,'temporary',
     1    l_temp2,k_temp2))
     2    call errquit('daxpyfile: MA problem',1)
        do i = 1,nblocks
          if (nodezero.and.util_print('daxpyfile',print_debug))
     1      write(LuOut,9210) (i-1)*blocksize,
     2                        min(size,i*blocksize)-(i-1)*blocksize
          k = 0
          do j = (i-1)*blocksize+1,min(size,i*blocksize)
            k = k + 1
            dbl_mb(k_temp1+k-1) = 0.0d0
          enddo
          do l = 1,n
            call get_block(d_in(l),dbl_mb(k_temp2),
     1        min(size,i*blocksize)-(i-1)*blocksize,(i-1)*blocksize)
            k = 0
            do j = (i-1)*blocksize+1,min(size,i*blocksize)
              k = k + 1
              dbl_mb(k_temp1+k-1) = dbl_mb(k_temp1+k-1) 
     1        + coefficients(l) * dbl_mb(k_temp2+k-1)
            enddo
          enddo
          call add_block(d_out,dbl_mb(k_temp1),
     1      min(size,i*blocksize)-(i-1)*blocksize,(i-1)*blocksize)
        enddo
        if (.not.ma_pop_stack(l_temp2))
     1    call errquit('createfile: MA problem',2)
        if (.not.ma_pop_stack(l_temp1))
     1    call errquit('createfile: MA problem',3)
      endif
      call ga_sync()
 9180 format(/,1x,'  input file handle:',i10,' size:',i10,
     1            ' coefficient:',f15.10)
 9190 format(/,1x,' output file handle:',i10,' size:',i10)
 9200 format(1x,'            nblocks:',i10,' size:',i10)
 9210 format(1x,'             offset:',i10,' size:',i10)
      return
      end
