      subroutine deletefile(d_sf)
c
c $Id: deletefile.F,v 1.4 2002-10-31 10:52:05 sohirata Exp $
c
c     Creates a shared file and returns the handle
c
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "util.fh"
#include "sf.fh"
#include "tce.fh"
      integer d_sf
      logical nodezero
      integer i,j
      logical used
      character*255 filename
      integer sf_unlink
      external sf_unlink
c
      nodezero = (ga_nodeid().eq.0)
      if (nodezero.and.util_print('files',print_debug))
     1  write(LuOut,9000) d_sf
      if (ioalg.eq.0) then
        call ga_sync()
        if (nfiles.eq.0)
     1    call errquit('deletefile: there is no file',nfiles)
        used = .false.
        do i = 1,nfiles
          if (filehandles(i).eq.d_sf) then
            filename = filenames(i)
            j = i
            used = .true.
          endif
        enddo
        if (.not.used)
     1    call errquit('deletefile: there is no such file',d_sf)
        if (j.lt.nfiles) then
          do i = j+1,nfiles
            filehandles(i-1) = filehandles(i)
            filenames(i-1) = filenames(i)
          enddo
        endif
        close(d_sf)
        call ga_sync()
        if (nodezero) then
          if (sf_unlink(filename).ne.0)
     1      call errquit('deletefile: cannot delete a file',0)
        endif
        nfiles = nfiles - 1
        if (nodezero.and.util_print('files',print_debug)) then
          do i = 1,nfiles
            write(LuOut,9010) filehandles(i),filenames(i)(1:40)
          enddo
        endif
        call ga_sync()
      else if (ioalg.eq.1) then
        if (sf_destroy(d_sf).ne.0)
     1    call errquit('deletefile: sf problem',0)
      else if (ioalg.eq.2) then
        if (.not.ga_destroy(d_sf))
     1    call errquit('deletefile: GA problem',0)
      endif
 9000 format(1x,'delete a file: file handle = ',i3)
 9010 format(1x,'file handle = ',i3,' file name = ',a)
      return
      end
