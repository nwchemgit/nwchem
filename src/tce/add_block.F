      subroutine add_block(d_file,array,size,offset)
c
c $Id$
c
c     This subroutine reads a section of a file and adds
c     the content of an array to it and write it back to
c     the original place in the file.  Note that two or more
c     different processes can read/add/write to the same 
c     section, the entire read/add/write must be protected
c     from interruption.
c
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "sf.fh"
#include "eaf.fh"
#include "util.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
      integer d_file
      integer d_f
      integer size
      integer offset
      double precision array(size)
      integer l_temp,k_temp
      integer request
      integer i,j
      logical used
      double precision element
      character*255 filename
      logical parallel
c
      cpusecs(4) = cpusecs(4) - util_cpusec()
      cpusecs(54) = cpusecs(54) - util_wallsec()
      parallel = (ga_nnodes().gt.1)
      if (util_print('add_block',print_debug)) then
         write(LuOut,9000) ga_nodeid(),d_file,size,offset
      endif
      if (offset .eq. -1)
     1  call errquit('add_block: illegal offset',0,UNKNOWN_ERR)
      if (ioalg.eq.0) then
        if (parallel) call ga_lock(0)
        used = .false.
        do i = 1,nfiles
          if (filehandles(i) .eq. d_file) then
            j = i
            used = .true.
          endif
        enddo
        if (.not.used) call errquit('add_block: invalid handle',d_file,
     1    UNKNOWN_ERR)
        if (parallel) open(d_file,file=filenames(j),
     1    access='direct',recl=bytes)
        do i = 1,size
          read(d_file,rec=offset+i) element
          write(d_file,rec=offset+i) array(i)+element
        enddo
        if (parallel) then
          close(d_file)
          call ga_unlock(0)
        endif
      else if (ioalg.eq.1) then
        if (.not.ma_push_get(mt_dbl,size,'temporary',
     1    l_temp,k_temp))
     2    call errquit('add_block: MA problem',0,MA_ERR)
        if (parallel) then
          used = .false.
          do i = 1,nfiles
            if (filehandles(i) .eq. d_file) then
              j = i
              used = .true.
            endif
          enddo
          filename = filenames(j)
          call ga_lock(0)
          if (eaf_open(filename, EAF_RW, d_f).ne.0)
     1      call errquit('add_block: eaf problem',2,DISK_ERR)
          if (eaf_read(d_f,dble(bytes)*dble(offset),
     1      dbl_mb(k_temp),bytes*size).ne.0)
     2      call errquit('add_block: eaf problem',2,DISK_ERR)
          do i = 1,size
            dbl_mb(k_temp+i-1) = dbl_mb(k_temp+i-1) + array(i)
          enddo
          if (eaf_write(d_f,dble(bytes)*dble(offset),
     1      dbl_mb(k_temp),bytes*size).ne.0)
     2      call errquit('add_block: eaf problem',2,DISK_ERR)
          if (eaf_close(d_f).ne.0)
     1      call errquit('add_block: eaf problem',0,DISK_ERR)
          call ga_unlock(0)
        else
          if (eaf_read(d_file,dble(bytes)*dble(offset),
     1      dbl_mb(k_temp),bytes*size).ne.0)
     2      call errquit('add_block: eaf problem',2,DISK_ERR)
          do i = 1,size
            dbl_mb(k_temp+i-1) = dbl_mb(k_temp+i-1) + array(i)
          enddo
          if (eaf_write(d_file,dble(bytes)*dble(offset),
     1      dbl_mb(k_temp),bytes*size).ne.0)
     2      call errquit('add_block: eaf problem',2,DISK_ERR)
        endif
        if (.not.ma_pop_stack(l_temp))
     1    call errquit('add_block: MA problem',0,MA_ERR)
      else if (ioalg.eq.2) then
        call ga_acc(d_file,offset+1,offset+size,1,1,
     1    array,1,1.0d0)
      else if (ioalg.eq.3) then
        if (.not.ma_push_get(mt_dbl,size,'temporary',
     1    l_temp,k_temp))
     2    call errquit('add_block: MA problem',0,MA_ERR)
        if (parallel) then
          call ga_lock(0)
          if (sf_open(d_file).ne.0)
     1      call errquit('add_block: sf problem',0,DISK_ERR)
        endif
        if (sf_read(d_file,dble(bytes)*dble(offset),
     1    dble(bytes)*dble(size),dbl_mb(k_temp),request).ne.0)
     2    call errquit('add_block: sf problem',1,DISK_ERR)
        if (sf_wait(request).ne.0)
     1    call errquit('add_block: sf problem',2,DISK_ERR)
        do i = 1,size
          dbl_mb(k_temp+i-1) = dbl_mb(k_temp+i-1) + array(i)
        enddo
        if (sf_write(d_file,dble(bytes)*dble(offset),
     1    dble(bytes)*dble(size),dbl_mb(k_temp),request).ne.0)
     2    call errquit('add_block: sf problem',3,DISK_ERR)
        if (sf_wait(request).ne.0)
     1    call errquit('add_block: sf problem',4,DISK_ERR)
        if (parallel) then
          if (sf_close(d_file).ne.0)
     1      call errquit('add_block: sf problem',5,DISK_ERR)
          call ga_unlock(0)
        endif
        if (.not.ma_pop_stack(l_temp))
     1    call errquit('add_block: MA problem',1,MA_ERR)
      else if (ioalg.eq.4) then
        if (.not.ma_push_get(mt_dbl,size,'temporary',
     1    l_temp,k_temp))
     2    call errquit('add_block: MA problem',0,MA_ERR)
        if (eaf_read(d_file,dble(bytes)*dble(offset),
     1    dbl_mb(k_temp),bytes*size).ne.0)
     2    call errquit('add_block: eaf problem',2,DISK_ERR)
        do i = 1,size
          dbl_mb(k_temp+i-1) = dbl_mb(k_temp+i-1) + array(i)
        enddo
        if (eaf_write(d_file,dble(bytes)*dble(offset),
     1    dbl_mb(k_temp),bytes*size).ne.0)
     2    call errquit('add_block: eaf problem',2,DISK_ERR)
        if (.not.ma_pop_stack(l_temp))
     1    call errquit('add_block: MA problem',0,MA_ERR)
      else if (ioalg.eq.5) then
        if (fileisga(d_file)) then
          call ga_acc(filehandles(d_file),offset+1,offset+size,1,1,
     1      array,1,1.0d0)
        else
          call errquit('add_block: read/write requests to DRA',0,
     1      UNKNOWN_ERR)
        endif
      else if (ioalg.eq.6) then
        if (fileisga(d_file)) then
          call ga_acc(filehandles(d_file),offset+1,offset+size,1,1,
     1      array,1,1.0d0)
        else
          if (.not.ma_push_get(mt_dbl,size,'temporary',
     1      l_temp,k_temp))
     2      call errquit('add_block: MA problem',0,MA_ERR)
          if (eaf_read(filehandles(d_file),
     1      dble(bytes)*dble(offset),
     2      dbl_mb(k_temp),bytes*size).ne.0)
     3      call errquit('add_block: eaf problem',2,DISK_ERR)
          do i = 1,size
            dbl_mb(k_temp+i-1) = dbl_mb(k_temp+i-1) + array(i)
          enddo
          if (eaf_write(filehandles(d_file),
     1      dble(bytes)*dble(offset),
     2      dbl_mb(k_temp),bytes*size).ne.0)
     3      call errquit('add_block: eaf problem',2,DISK_ERR)
          if (.not.ma_pop_stack(l_temp))
     1      call errquit('add_block: MA problem',0,MA_ERR)
        endif
      endif
 9000 format(1x,'node',i3,' add_block request to file:',i10,
     1  ' size:',i10,' offset:',i10)
      cpusecs(4) = cpusecs(4) + util_cpusec()
      cpusecs(54) = cpusecs(54) + util_wallsec()
      return
      end












      subroutine add_block_nb(d_file,array,size,offset,nbh)
c
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "util.fh"
#include "tce.fh"
#include "tce_main.fh"
c
      integer d_file
      integer size
      integer offset
      integer nbh
      double precision array(size)
      logical parallel
c
      cpusecs(4) = cpusecs(4) - util_cpusec()
      cpusecs(54) = cpusecs(54) - util_wallsec()
c
      parallel = (ga_nnodes().gt.1)
c
      if (util_print('add_block',print_debug)) then
         write(LuOut,9000) ga_nodeid(),d_file,size,offset
      endif
c
      if (offset .eq. -1)
     1  call errquit('add_block_nb: illegal offset',0,UNKNOWN_ERR)
c
      if (ioalg.eq.2) then
        call ga_nbacc(d_file,offset+1,offset+size,1,1,array,1,1.0d0,nbh)
      else
        call errquit('add_block_nb: only for GA!!!',911,MA_ERR)
      endif
c
      cpusecs(4) = cpusecs(4) + util_cpusec()
      cpusecs(54) = cpusecs(54) + util_wallsec()
      return
c
 9000 format(1x,'node',i3,' add_block_nb request to file:',i10,
     1  ' size:',i10,' offset:',i10)
c
      end
