      subroutine put_block(d_file,array,size,offset)
c
c $Id: put_block.F,v 1.10 2002-12-16 21:13:38 sohirata Exp $
c
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "sf.fh"
#include "util.fh"
#include "stdio.fh"
#include "tce.fh"
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
     1  call errquit('put_block: illegal offset',0)
      if (ioalg.eq.0) then
        used = .false.
        do i = 1,nfiles
          if (filehandles(i) .eq. d_file) then
            j = i
            used = .true.
          endif
        enddo
        if (.not.used) call errquit('put_block: invalid handle',d_file)
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
     1      call errquit('put_block: invalid handle',d_file)
          call ga_lock(0)
          d_f = f_open(filename)
          if (f_write(d_f,bytes*offset,bytes*size,array).ne.0)
     1      call errquit('put_block: file problem',0)
          if (f_close(d_f).ne.0)
     1      call errquit('put_block: file problem',1)
          call ga_unlock(0)
        else
          if (f_write(d_file,bytes*offset,bytes*size,array).ne.0)
     1      call errquit('put_block: file problem',0)
        endif
      else if (ioalg.eq.2) then
        call ga_put(d_file,offset+1,offset+size,1,1,
     1    array,1)
      else if (ioalg.eq.3) then
        if (parallel) then
          call ga_lock(0)
          if (sf_open(d_file).ne.0)
     1      call errquit('put_block: sf problem',0)
        endif
        if (sf_write(d_file,dfloat(bytes)*dfloat(offset),
     1    dfloat(bytes)*dfloat(size),array,request).ne.0)
     2    call errquit('put_block: sf problem',1)
        if (sf_wait(request).ne.0)
     1    call errquit('put_block: sf problem',2)
        if (parallel) then
          if (sf_close(d_file).ne.0)
     1      call errquit('put_block: sf problem',3)
          call ga_unlock(0)
        endif
      else if (ioalg.eq.4) then
        if (f_write(d_file,bytes*offset,bytes*size,array).ne.0)
     1    call errquit('put_block: file problem',0)
      endif
 9000 format(1x,'node',i3,' put_block request to file:',i10,
     1  ' size:',i10,' offset:',i10)
      cpusecs(3) = cpusecs(3) + util_cpusec()
      cpusecs(53) = cpusecs(53) + util_wallsec()
      return
      end
