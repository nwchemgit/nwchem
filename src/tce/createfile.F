      subroutine createfile(filename,d_sf,size)
c
c $Id: createfile.F,v 1.6 2002-10-31 10:52:05 sohirata Exp $
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
      character*(*) filename
      integer d_sf
      integer size
      integer i
      logical nodezero
      logical used
c
      nodezero = (ga_nodeid().eq.0)
      if (nodezero.and.util_print('files',print_debug))
     1  write(LuOut,9000) size,filename(1:40)
      if (ioalg.eq.0) then
        call ga_sync()
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
        call errquit('createfile: file system problem',0)
  100   continue
        if (nodezero.and.util_print('files',print_debug)) then
          do i = 1,nfiles
            write(LuOut,9010) filehandles(i),filenames(i)(1:40)
          enddo
        endif
        call ga_lock(0)
        open(d_sf,file=filenames(nfiles),form='unformatted',
     1    access='direct',recl=bytes)
        close(d_sf)
        call ga_unlock(0)
        call ga_sync()
        call tce_zero(d_sf,size)
      else if (ioalg.eq.1) then
        if (sf_create(filename,dfloat(bytes*size),dfloat(bytes*size),
     1    1.0d4,d_sf).ne.0)
     2    call errquit('createfile: sf problem',0)
        call tce_zero(d_sf,size)
      else if (ioalg.eq.2) then
        if (.not.ga_create(mt_dbl,size,1,filename,
     1    -1,1,d_sf)) call errquit
     2    ('createfile: GA problem',0)
        call tce_zero(d_sf,size)
      endif
 9000 format(1x,'create a file: size = ',i10,' file name = ',a)
 9010 format(1x,'file handle = ',i3,' file name = ',a)
      return
      end
