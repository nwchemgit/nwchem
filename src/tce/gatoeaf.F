      subroutine gatoeaf(d_sf)
c
c $Id$
c
c     Unload data from GA to EAF
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
      integer i
      integer size
      logical nodezero
      integer nblocks
      integer blocksize
      integer l_temp, k_temp
c
      nodezero = (ga_nodeid().eq.0)
      if (ioalg.eq.6) then
        if (.not.fileisga(d_sf))
     1    call errquit('gatoeaf: no GA found',0,UNKNOWN_ERR)
        call ga_sync()
        size = filesizes(d_sf)
        if (eaf_open(filenames(d_sf), EAF_RW, d_temp).ne.0)
     1    call errquit('gatoeaf: eaf problem',0,DISK_ERR)
        nblocks = size/buffer_size
        if (nblocks*buffer_size.lt.size) nblocks = nblocks + 1
        blocksize = size/nblocks
        if (nblocks*blocksize.lt.size) blocksize = blocksize + 1
        if (.not.ma_push_get(mt_dbl,blocksize,'temporary',
     1    l_temp,k_temp))
     2    call errquit('gatoeaf: MA problem',0,MA_ERR)
        do i = 1,nblocks
          call ga_get(filehandles(d_sf),(i-1)*blocksize+1,
     1      min(size,i*blocksize),1,1,dbl_mb(k_temp),1)
          if (eaf_write(d_temp,dble(bytes)*dble((i-1)*blocksize),
     1      dbl_mb(k_temp),
     2      bytes*(min(size,i*blocksize)-(i-1)*blocksize)).ne.0)
     3      call errquit('gatoeaf: eaf problem',7,DISK_ERR)
        enddo
        if (.not.ma_pop_stack(l_temp))
     1    call errquit('gatoeaf: MA problem',1,MA_ERR)
        if (.not.ga_destroy(filehandles(d_sf)))
     1    call errquit('gatoeaf: GA problem',0,GA_ERR)
        filehandles(d_sf) = d_temp
        fileisga(d_sf) = .false.
        call ga_sync()
      endif
      return
      end
