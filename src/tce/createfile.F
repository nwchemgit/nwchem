      subroutine createfile(filename,d_sf,size)
c
c $Id: createfile.F,v 1.15 2003-09-16 21:50:11 sohirata Exp $
c
c     Creates a shared file and returns the handle
c
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "util.fh"
#include "sf.fh"
#include "eaf.fh"
#include "dra.fh"
#include "tce.fh"
      character*(*) filename
      integer d_sf
      integer d_temp
      integer size
      integer i
      logical nodezero
      logical parallel
      logical used
      integer request
c
      cpusecs(1) = cpusecs(1) - util_cpusec()
      cpusecs(51) = cpusecs(51) - util_wallsec()
      nodezero = (ga_nodeid().eq.0)
      parallel = (ga_nnodes().gt.1)
      if (ioalg.eq.0) then
        if (nodezero.and.util_print('files',print_debug))
     1    write(LuOut,9000) size,filename(1:40)
        if (nfiles.eq.max_files)
     1    call errquit('createfile: exceeded max_files',nfiles)
        do d_sf = 10,max_files+9
          used = .false.
          do i = 1,nfiles
            if (filehandles(i).eq.d_sf) used = .true.
          enddo
          if (.not.used) then
            nfiles = nfiles + 1
            filehandles(nfiles) = d_sf
            filenames(nfiles) = filename
            goto 100
          endif
        enddo
        call errquit('createfile: file problem',0)
  100   continue
        if (nodezero.and.util_print('files',print_debug)) then
          do i = 1,nfiles
            write(LuOut,9010) filehandles(i),filenames(i)(1:40)
          enddo
        endif
        if (parallel) then
          call ga_sync()
          call ga_lock(0)
        endif
        open(d_sf,file=filenames(nfiles),form='unformatted',
     1    access='direct',recl=bytes)
        if (parallel) then
          close(d_sf)
          call ga_unlock(0)
          call ga_sync()
        endif
        call tce_zero(d_sf,size)
      else if (ioalg.eq.1) then
        if (nodezero.and.util_print('files',print_debug))
     1    write(LuOut,9000) size,filename(1:40)
        if (parallel) then
          call ga_sync()
          if (nfiles.eq.max_files)
     1      call errquit('createfile: exceeded max_files',nfiles)
          do d_sf = 10,max_files+9
            used = .false.
            do i = 1,nfiles
              if (filehandles(i).eq.d_sf) used = .true.
            enddo
            if (.not.used) then
              nfiles = nfiles + 1
              filehandles(nfiles) = d_sf
              filenames(nfiles) = filename
              goto 110
            endif
          enddo
          call errquit('createfile: file problem',0)
  110     continue
          if (eaf_open(filename, EAF_RW, d_temp).ne.0)
     1      call errquit('createfile: eaf problem',0)
          if (eaf_close(d_temp).ne.0)
     1      call errquit('createfile: eaf problem',1)
          if (nodezero.and.util_print('files',print_debug)) then
            do i = 1,nfiles
              write(LuOut,9010) filehandles(i),filenames(i)(1:40)
            enddo
          endif
        else
          if (nfiles.eq.max_files)
     1      call errquit('createfile: exceeded max_files',nfiles)
          if (eaf_open(filename, EAF_RW, d_sf).ne.0)
     1      call errquit('createfile: eaf problem',0)
          nfiles = nfiles + 1
          filehandles(nfiles) = d_sf
          filenames(nfiles) = filename
        endif
        call tce_zero(d_sf,size)
      else if (ioalg.eq.2) then
        if (.not.ga_create(mt_dbl,size,1,filename,
     1    -1,1,d_sf)) then
          write(LuOut,*) ' available GA memory ',
     1    ga_memory_avail(),' bytes'
          call errquit ('createfile: failed ga_create size=',size)
        endif
        if (nodezero.and.util_print('files',print_debug))
     1    write(LuOut,9010) d_sf,filename(1:40)
        call tce_zero(d_sf,size)
      else if (ioalg.eq.3) then
        if (parallel) call ga_sync()
        if (sf_create(filename,dfloat(bytes)*dfloat(size),
     1      dfloat(bytes)*dfloat(size),1.0d4,d_sf).ne.0)
     2    call errquit('createfile: sf problem',0)
        if (parallel) then
          if (sf_close(d_sf).ne.0)
     1      call errquit('createfile: sf problem',1)
          call ga_sync()
        endif
        if (nodezero.and.util_print('files',print_debug))
     1    write(LuOut,9010) d_sf,filename(1:40)
        call tce_zero(d_sf,size)
      else if (ioalg.eq.4) then
        if (nfiles.eq.max_files)
     1    call errquit('createfile: exceeded max_files',nfiles)
        if (eaf_open(filename, EAF_RW, d_sf).ne.0)
     1    call errquit('createfile: eaf problem',0)
        nfiles = nfiles + 1
        filehandles(nfiles) = d_sf
        filenames(nfiles) = filename
        if (util_print('files',print_debug))
     1    write(LuOut,9010) d_sf,filename(1:40)
        call tce_zero(d_sf,size)
      else if (ioalg.eq.5) then
        if (nfiles.eq.max_files)
     1    call errquit('createfile: exceeded max_files',nfiles)
        if (.not.ga_create(mt_dbl,size,1,filename,
     1    -1,1,d_temp))
     2    call errquit ('createfile: ga problem',1)
        call ga_zero(d_temp)
        nfiles = nfiles + 1
        filehandles(nfiles) = d_temp
        filenames(nfiles) = filename
        filesizes(nfiles) = size
        fileisga(nfiles) = .true.
        d_sf = nfiles
      else if (ioalg.eq.6) then
        if (nfiles.eq.max_files)
     1    call errquit('createfile: exceeded max_files',nfiles)
        if (.not.ga_create(mt_dbl,size,1,filename,
     1    -1,1,d_temp))
     2    call errquit ('createfile: ga problem',1)
        nfiles = nfiles + 1
        filehandles(nfiles) = d_temp
        filenames(nfiles) = filename
        filesizes(nfiles) = size
        fileisga(nfiles) = .true.
        d_sf = nfiles
        call tce_zero(d_sf,size)
        if (util_print('files',print_debug))
     1    call filesummary
      endif
      cpusecs(1) = cpusecs(1) + util_cpusec()
      cpusecs(51) = cpusecs(51) + util_wallsec()
 9000 format(1x,'create a file: size = ',i10,' file name = ',a)
 9010 format(1x,'file handle = ',i9,' file name = ',a)
 9020 format(1x,'node',i3,' file handle = ',i9,' file name = ',a)
 9030 format(1x,'create a file: node',i3,' size = ',i10,
     1          ' file name = ',a)
      return
      end
