      subroutine deletefile(d_sf)
c
c $Id$
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
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
      integer d_sf
      integer d_temp
      logical nodezero
      logical parallel
      integer i,j
      logical used
      character*255 filename
c
      cpusecs(2) = cpusecs(2) - util_cpusec()
      cpusecs(52) = cpusecs(52) - util_wallsec()
      nodezero = (ga_nodeid().eq.0)
      parallel = (ga_nnodes().gt.1)
      if (nodezero.and.util_print('files',print_debug))
     1  write(LuOut,9000) d_sf
      if (ioalg.eq.0) then
        if (parallel) call ga_sync()
        if (nfiles.eq.0)
     1    call errquit('deletefile: there is no file',nfiles,
     2    UNKNOWN_ERR)
        used = .false.
        do i = 1,nfiles
          if (filehandles(i).eq.d_sf) then
            filename = filenames(i)
            j = i
            used = .true.
          endif
        enddo
        if (.not.used)
     1    call errquit('deletefile: there is no such file',d_sf,
     2    UNKNOWN_ERR)
        if (j.lt.nfiles) then
          do i = j+1,nfiles
            filehandles(i-1) = filehandles(i)
            filenames(i-1) = filenames(i)
          enddo
        endif
        close(d_sf)
        if (parallel) call ga_sync()
        if (nodezero) then
          call util_file_unlink(filename)
        endif
        nfiles = nfiles - 1
        if (nodezero.and.util_print('files',print_debug)) then
          do i = 1,nfiles
            write(LuOut,9010) filehandles(i),filenames(i)(1:40)
          enddo
        endif
        if (parallel) call ga_sync()
      else if (ioalg.eq.1) then
        if (nfiles.eq.0)
     1    call errquit('deletefile: there is no file',nfiles,
     2    UNKNOWN_ERR)
        used = .false.
        do i = 1,nfiles
          if (filehandles(i).eq.d_sf) then
            filename = filenames(i)
            j = i
            used = .true.
          endif
        enddo
        if (.not.used)
     1    call errquit('deletefile: there is no such file',d_sf,
     2    UNKNOWN_ERR)
        if (j.lt.nfiles) then
          do i = j+1,nfiles
            filehandles(i-1) = filehandles(i)
            filenames(i-1) = filenames(i)
          enddo
        endif
        if (parallel) then
          call ga_sync()
          if (nodezero) then
            call util_file_unlink(filename)
          endif
          call ga_sync()
        else
          if (eaf_close(d_sf).ne.0)
     1      call errquit('deletefile: eaf problem',0,DISK_ERR)
          call util_file_unlink(filename)
        endif
        nfiles = nfiles - 1
        if (nodezero.and.util_print('files',print_debug)) then
          do i = 1,nfiles
            write(LuOut,9010) filehandles(i),filenames(i)(1:40)
          enddo
        endif
      else if (ioalg.eq.2) then
        if (.not.ga_destroy(d_sf))
     1    call errquit('deletefile: GA problem',0,GA_ERR)
      else if (ioalg.eq.3) then
        if (parallel) then
          call ga_sync()
          if (sf_open(d_sf).ne.0)
     1      call errquit('deletefile: sf problem',0,DISK_ERR)
        endif
        if (sf_destroy(d_sf).ne.0)
     1    call errquit('deletefile: sf problem',1,DISK_ERR)
        if (parallel) call ga_sync()
      else if (ioalg.eq.4) then
        if (nfiles.eq.0)
     1    call errquit('deletefile: there is no file',nfiles,
     2    UNKNOWN_ERR)
        used = .false.
        do i = 1,nfiles
          if (filehandles(i).eq.d_sf) then
            filename = filenames(i)
            j = i
            used = .true.
          endif
        enddo
        if (.not.used)
     1    call errquit('deletefile: there is no such file',d_sf,
     2    UNKNOWN_ERR)
        if (j.lt.nfiles) then
          do i = j+1,nfiles
            filehandles(i-1) = filehandles(i)
            filenames(i-1) = filenames(i)
          enddo
        endif
        if (eaf_close(d_sf).ne.0)
     1    call errquit('deletefile: eaf problem',0,DISK_ERR)
        call util_file_unlink(filename)
        nfiles = nfiles - 1
        if (nodezero.and.util_print('files',print_debug)) then
          do i = 1,nfiles
            write(LuOut,9010) filehandles(i),filenames(i)(1:40)
          enddo
        endif
      else if (ioalg.eq.5) then
        if (fileisga(d_sf)) then
          if (.not.ga_destroy(filehandles(d_sf)))
     1      call errquit('deletefile: GA problem',0,GA_ERR)
        else
          if (dra_open(filenames(d_sf),DRA_R,d_temp).ne.0)
     1      call errquit('deletefile: dra problem',0,DISK_ERR)
          if (dra_delete(d_temp).ne.0)
     1      call errquit('deletefile: dra problem',1,DISK_ERR)
        endif
        filehandles(d_sf) = 2**30
        if (d_sf.eq.nfiles) nfiles = nfiles - 1
        if (nodezero.and.util_print('files',print_debug))
     1    call filesummary
      else if (ioalg.eq.6) then
        if (fileisga(d_sf)) then
          if (.not.ga_destroy(filehandles(d_sf))) then
            write(LuOut,9010) filehandles(d_sf),filenames(d_sf)
            call errquit('deletefile: GA problem',0,GA_ERR)
          endif
        else
          if (eaf_close(filehandles(d_sf)).ne.0) then
            write(LuOut,9010) filehandles(d_sf),filenames(d_sf)
            call errquit('deletefile: eaf problem',0,DISK_ERR)
          endif
          call util_file_unlink(filenames(d_sf))
        endif
        filehandles(d_sf) = 2**30
        if (d_sf.eq.nfiles) nfiles = nfiles - 1
      endif
      cpusecs(2) = cpusecs(2) + util_cpusec()
      cpusecs(52) = cpusecs(52) + util_wallsec()
 9000 format(1x,'delete a file: file handle = ',i9)
 9010 format(1x,'file handle = ',i9,' file name = ',a)
      return
      end



















      subroutine deletefile_sf(d_sf)
c
c $Id$
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
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
      integer d_sf
      integer d_temp
      logical nodezero
      logical parallel
      integer i,j
      logical used
      character*255 filename
c
      cpusecs(2) = cpusecs(2) - util_cpusec()
      cpusecs(52) = cpusecs(52) - util_wallsec()
      nodezero = (ga_nodeid().eq.0)
      parallel = (ga_nnodes().gt.1)
      if (nodezero.and.util_print('files',print_debug))
     1  write(LuOut,9000) d_sf
c
c      else if (ioalg.eq.3) then
c
        if (parallel) then
          call ga_sync()
          if (sf_open(d_sf).ne.0)
     1      call errquit('deletefile: sf problem',0,DISK_ERR)
        endif
        if (sf_destroy(d_sf).ne.0)
     1    call errquit('deletefile: sf problem',1,DISK_ERR)
        if (parallel) call ga_sync()
c
c      endif
c
      cpusecs(2) = cpusecs(2) + util_cpusec()
      cpusecs(52) = cpusecs(52) + util_wallsec()
 9000 format(1x,'delete a file: file handle = ',i9)
 9010 format(1x,'file handle = ',i9,' file name = ',a)
      return
      end
