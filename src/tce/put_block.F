      subroutine put_block(d_file,array,size,offset)
c
c $Id$
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
      integer request
      integer i,j
      logical used
      character*255 filename
      logical parallel
c
      cpusecs(3) = cpusecs(3) - util_cpusec()
      cpusecs(53) = cpusecs(53) - util_wallsec()
      parallel = (ga_nnodes().gt.1)
      if (util_print('put_block',print_debug)) then
         write(LuOut,9000) ga_nodeid(),d_file,size,offset
      endif
      if (offset .eq. -1)
     1  call errquit('put_block: illegal offset',0,UNKNOWN_ERR)
      if (ioalg.eq.0) then
        used = .false.
        do i = 1,nfiles
          if (filehandles(i) .eq. d_file) then
            j = i
            used = .true.
          endif
        enddo
        if (.not.used) call errquit('put_block: invalid handle',d_file,
     1  UNKNOWN_ERR)
        if (parallel) then
          call ga_lock(0)
          open(d_file,file=filenames(j),access='direct',recl=bytes)
        endif
        do i = 1,size
          write(d_file,rec=offset+i) array(i)
        enddo
        if (parallel) then
          close(d_file)
          call ga_unlock(0)
        endif
      else if (ioalg.eq.1) then
        if (parallel) then
          used = .false.
          do i = 1,nfiles
            if (filehandles(i) .eq. d_file) then
              j = i
              used = .true.
            endif
          enddo
          filename = filenames(j)
          if (.not.used)
     1      call errquit('put_block: invalid handle',d_file,UNKNOWN_ERR)
          call ga_lock(0)
          if (eaf_open(filename, EAF_RW, d_f).ne.0)
     1      call errquit('put_block: eaf problem',1,DISK_ERR)
          if (eaf_write(d_f,dble(bytes)*dble(offset),
     1      array,bytes*size).ne.0)
     2      call errquit('put_block: eaf problem',2,DISK_ERR)
          if (eaf_close(d_f).ne.0)
     1      call errquit('put_block: eaf problem',4,DISK_ERR)
          call ga_unlock(0)
        else
          if (eaf_write(d_file,dble(bytes)*dble(offset),
     1      array,bytes*size).ne.0)
     2      call errquit('put_block: eaf problem',5,DISK_ERR)
        endif
      else if (ioalg.eq.2) then
        call ga_put(d_file,offset+1,offset+size,1,1,
     1    array,size)
      else if (ioalg.eq.3) then
        if (parallel) then
          call ga_lock(0)
          if (sf_open(d_file).ne.0)
     1      call errquit('put_block: sf problem',0,DISK_ERR)
        endif
        if (sf_write(d_file,dble(bytes)*dble(offset),
     1    dble(bytes)*dble(size),array,request).ne.0)
     2    call errquit('put_block: sf problem',1,DISK_ERR)
        if (sf_wait(request).ne.0)
     1    call errquit('put_block: sf problem',2,DISK_ERR)
        if (parallel) then
          if (sf_close(d_file).ne.0)
     1      call errquit('put_block: sf problem',3,DISK_ERR)
          call ga_unlock(0)
        endif
      else if (ioalg.eq.4) then
        if (eaf_write(d_file,dble(bytes)*dble(offset),
     1    array,bytes*size).ne.0)
     2    call errquit('put_block: eaf problem',2,DISK_ERR)
      else if (ioalg.eq.5) then
        if (fileisga(d_file)) then
          call ga_put(filehandles(d_file),offset+1,offset+size,1,1,
     1      array,size)
        else
          call errquit('put_block: a write request to DRA',1,
     1      UNKNOWN_ERR)
        endif
      else if (ioalg.eq.6) then
        if (fileisga(d_file)) then
          call ga_put(filehandles(d_file),offset+1,offset+size,1,1,
     1      array,size)
        else
          if (eaf_write(filehandles(d_file),
     1      dble(bytes)*dble(offset),array,bytes*size).ne.0)
     2      call errquit('put_block: eaf problem',2,DISK_ERR)
        endif
      endif
 9000 format(1x,'node',i3,' put_block request to file:',i10,
     1  ' size:',i10,' offset:',i10)
      cpusecs(3) = cpusecs(3) + util_cpusec()
      cpusecs(53) = cpusecs(53) + util_wallsec()
      return
      end




















      subroutine put_block_sf(d_file,array,size,offset)
c
c $Id$
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
      integer request
      integer i,j
      logical used
      character*255 filename
      logical parallel
c
      cpusecs(3) = cpusecs(3) - util_cpusec()
      cpusecs(53) = cpusecs(53) - util_wallsec()
      parallel = (ga_nnodes().gt.1)
      if (util_print('put_block',print_debug)) then
         write(LuOut,9000) ga_nodeid(),d_file,size,offset
      endif
      if (offset .eq. -1)
     1  call errquit('put_block: illegal offset',0,UNKNOWN_ERR)
c
c      else if (ioalg.eq.3) then
c
        if (parallel) then
          call ga_lock(0)
          if (sf_open(d_file).ne.0)
     1      call errquit('put_block: sf problem',0,DISK_ERR)
        endif
        if (sf_write(d_file,dble(bytes)*dble(offset),
     1    dble(bytes)*dble(size),array,request).ne.0)
     2    call errquit('put_block: sf problem',1,DISK_ERR)
        if (sf_wait(request).ne.0)
     1    call errquit('put_block: sf problem',2,DISK_ERR)
        if (parallel) then
          if (sf_close(d_file).ne.0)
     1      call errquit('put_block: sf problem',3,DISK_ERR)
          call ga_unlock(0)
        endif
c
c      endif
c
 9000 format(1x,'node',i3,' put_block request to file:',i10,
     1  ' size:',i10,' offset:',i10)
      cpusecs(3) = cpusecs(3) + util_cpusec()
      cpusecs(53) = cpusecs(53) + util_wallsec()
      return
      end
