      subroutine reconcilefile(d_sf,size)
c
c $Id: reconcilefile.F,v 1.5 2002-12-17 00:03:10 sohirata Exp $
c
c     Accumulate local distributed files and update all files
c     with the one accumulated file
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
      integer g_temp
      integer l_temp,k_temp
      integer i
      logical nodezero
      logical parallel
c 
      nodezero = (ga_nodeid().eq.0)
      parallel = (ga_nnodes().gt.1)
      if (.not. parallel) return
      if (ioalg.ne.4) then
        call ga_sync()
        return
      endif
      cpusecs(6) = cpusecs(6) - util_cpusec()
      cpusecs(56) = cpusecs(56) - util_wallsec()
      if (util_print('files',print_debug))
     1  write(LuOut,9190) ga_nodeid(),d_sf,size
      nblocks = size/buffer_size + 1
      blocksize = size/nblocks + 1
      if (.not.ga_create(mt_dbl,blocksize,1,'temporary',
     1  -1,1,g_temp))
     2  call errquit('reconcilefile: GA problem',0)
      if (.not.ma_push_get(mt_dbl,blocksize,'temporary',
     1  l_temp,k_temp))
     2  call errquit('reconcilefile: MA problem',0)
      do i = 1,nblocks
        call ga_zero(g_temp)
        call get_block(d_sf,dbl_mb(k_temp),
     1    min(size,i*blocksize)-(i-1)*blocksize,(i-1)*blocksize)
        call ga_acc(g_temp,1,blocksize,1,1,dbl_mb(k_temp),1,1.0d0)
        call ga_sync()
        call ga_get(g_temp,1,blocksize,1,1,dbl_mb(k_temp),1)
        call put_block(d_sf,dbl_mb(k_temp),
     1    min(size,i*blocksize)-(i-1)*blocksize,(i-1)*blocksize)
      enddo
      call ga_sync()
      if (.not.ga_destroy(g_temp))
     1  call errquit('reconcilefile: GA problem',0)
      if (.not.ma_pop_stack(l_temp))
     1  call errquit('reconcilefile: MA problem',1)
 9190 format(/,1x,'Reconcile file: node',i3,' handle:',i10,
     1            ' size:',i10)
 9200 format(1x,'        offset:',i10,' size:',i10)
      cpusecs(6) = cpusecs(6) + util_cpusec()
      cpusecs(56) = cpusecs(56) + util_wallsec()
      return
      end
