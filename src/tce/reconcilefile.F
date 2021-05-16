      subroutine reconcilefile(d_sf,size)
c
c $Id$
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
#include "eaf.fh"
#include "dra.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
      integer d_sf
      integer d_temp
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
      cpusecs(6) = cpusecs(6) - util_cpusec()
      cpusecs(56) = cpusecs(56) - util_wallsec()
      if (ioalg.eq.0) then
        if (parallel) call ga_sync()
      else if (ioalg.eq.1) then
        if (parallel) call ga_sync()
      else if (ioalg.eq.2) then
        if (parallel) call ga_sync()
        if (1.eq.ga_is_mirrored(d_sf)) then
          if (util_print('reconcilefile',print_debug)) write(LuOut,*) 
     1        'accumulating mirrored array across all nodes'
          call ga_merge_mirrored(d_sf)
        endif
      else if (ioalg.eq.3) then
        if (parallel) call ga_sync()
      else if ((ioalg.eq.4).or.
     1        ((ioalg.eq.6).and.(.not.fileisga(d_sf)))) then
        if (parallel) then
          if (util_print('reconcilefile',print_debug))
     1      write(LuOut,9190) ga_nodeid(),d_sf,size
          if (nodezero.and.(ioalg.eq.6)) write(LuOut,9230)
          nblocks = size/buffer_size
          if (nblocks*buffer_size.lt.size) nblocks = nblocks + 1
          blocksize = size/nblocks
          if (nblocks*blocksize.lt.size) blocksize = blocksize + 1
          if (.not.ga_create(mt_dbl,blocksize,1,'temporary',
     1      -1,1,g_temp))
     2      call errquit('reconcilefile: GA problem',0,GA_ERR)
          if (.not.ma_push_get(mt_dbl,blocksize,'temporary',
     1      l_temp,k_temp))
     2      call errquit('reconcilefile: MA problem',0,MA_ERR)
          do i = 1,nblocks
            call ga_zero(g_temp)
            call get_block(d_sf,dbl_mb(k_temp),
     1        min(size,i*blocksize)-(i-1)*blocksize,(i-1)*blocksize)
            call ga_acc(g_temp,1,blocksize,1,1,dbl_mb(k_temp),1,1.0d0)
            call ga_sync()
            call ga_get(g_temp,1,blocksize,1,1,dbl_mb(k_temp),1)
            call put_block(d_sf,dbl_mb(k_temp),
     1        min(size,i*blocksize)-(i-1)*blocksize,(i-1)*blocksize)
          enddo
          call ga_sync()
          if (.not.ga_destroy(g_temp))
     1      call errquit('reconcilefile: GA problem',0,GA_ERR)
          if (.not.ma_pop_stack(l_temp))
     1      call errquit('reconcilefile: MA problem',1,MA_ERR)
          if (nodezero.and.util_print('reconcilefile',print_debug))
     1      call filesummary
        endif
      else if (ioalg.eq.5) then
        if (parallel) call ga_sync()
      else if ((ioalg.eq.6).and.fileisga(d_sf)) then
        call ga_sync()
        if (eaf_open(filenames(d_sf), EAF_RW, d_temp).ne.0)
     1    call errquit('reconcilefile: eaf problem',0,DISK_ERR)
        if (nodezero.and.util_print('reconcilefile',print_debug))
     1    write(LuOut,9210) d_sf,filehandles(d_sf),d_temp
        nblocks = size/buffer_size
        if (nblocks*buffer_size.lt.size) nblocks = nblocks + 1
        blocksize = size/nblocks
        if (nblocks*blocksize.lt.size) blocksize = blocksize + 1
        if (.not.ma_push_get(mt_dbl,blocksize,'temporary',
     1    l_temp,k_temp))
     2    call errquit('reconcilefile: MA problem',0,MA_ERR)
        do i = 1,nblocks
          call ga_get(filehandles(d_sf),(i-1)*blocksize+1,
     1      min(size,i*blocksize),1,1,dbl_mb(k_temp),1)
          if (eaf_write(d_temp,dble(bytes)*dble((i-1)*blocksize),
     1      dbl_mb(k_temp),
     2      bytes*(min(size,i*blocksize)-(i-1)*blocksize)).ne.0)
     3      call errquit('reconcilefile: eaf problem',7,DISK_ERR)
        enddo
        if (.not.ma_pop_stack(l_temp))
     1    call errquit('reconcilefile: MA problem',1,MA_ERR)
        if (.not.ga_destroy(filehandles(d_sf)))
     1    call errquit('reconcilefile: GA problem',0,GA_ERR)
        filehandles(d_sf) = d_temp
        fileisga(d_sf) = .false.
        if (nodezero.and.util_print('reconcilefile',print_debug))
     1    call filesummary
        call ga_sync()
      endif
      cpusecs(6) = cpusecs(6) + util_cpusec()
      cpusecs(56) = cpusecs(56) + util_wallsec()
 9190 format(/,1x,'Reconcile file: node',i3,' handle:',i10,
     1            ' size:',i10)
 9200 format(1x,'        offset:',i10,' size:',i10)
 9210 format(/,1x,'Reconcile file: logical handle:',i10,
     1            ' GA handle:',i10,' EAF handle:',i10)
 9220 format(/,1x,i3,' reconcile file: logical handle:',i10,
     1            ' EAF handle:',i10)
 9230 format(/,1x,'Warning: EAF reconciliation in GA/EAF algorithm')
      return
 9240 format(/,1x,'Reconcile file: logical handle:',i10,
     1            ' GA handle:',i10,' DRA handle:',i10)
      end



      subroutine filesummary
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "util.fh"
#include "sf.fh"
#include "eaf.fh"
#include "tce.fh"
#include "tce_main.fh"
      integer i
      logical nodezero
      nodezero = (ga_nodeid().eq.0)
      if (nodezero) then
        do i = 1,nfiles
          if (ioalg.eq.5) then
            if (fileisga(i)) write(LuOut,9000) i,filehandles(i),
     1        filenames(i)(1:40)
            if (.not.fileisga(i)) write(LuOut,9010) i,filehandles(i),
     1        filenames(i)(1:40)
          else if (ioalg.eq.6) then
            if (fileisga(i)) write(LuOut,9000) i,filehandles(i),
     1        filenames(i)(1:40)
            if (.not.fileisga(i)) write(LuOut,9020) i,filehandles(i),
     1        filenames(i)(1:40)
          endif
        enddo
      endif
 9000 format(1x,i3,' GA  file handle = ',i9,' file name = ',a)
 9010 format(1x,i3,' DRA file handle = ',i9,' file name = ',a)
 9020 format(1x,i3,' EAF file handle = ',i9,' file name = ',a)
      return
      end
