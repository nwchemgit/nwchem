      subroutine createfile(filename,d_sf,size)
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
      character*(*) filename
      integer d_sf
      integer d_temp
      integer size
      integer i
      logical nodezero
      logical parallel
      logical used
c
      nodezero = (ga_nodeid().eq.0)
      parallel = (ga_nnodes().gt.1)
      cpusecs(1) = cpusecs(1) - util_cpusec()
      cpusecs(51) = cpusecs(51) - util_wallsec()
      if (size.eq.0) then
        if (nodezero.and.util_print('files',print_high)) then
          write(LuOut,*) '*************** Warning ***************'
          write(LuOut,*) ' Create file size is zero. Calculation'
          write(LuOut,*) ' will continue by increasing the size.'
          write(LuOut,*) ' Use of a larger basis set is advised.'
        endif
        size = 1
      endif
      if (ioalg.eq.0) then
        if (nodezero.and.util_print('files',print_debug))
     1    write(LuOut,9000) size,filename(1:40)
        if (nfiles.eq.max_files)
     1    call errquit('createfile: exceeded max_files',nfiles,
     2    UNKNOWN_ERR)
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
        call errquit('createfile: file problem',0,UNKNOWN_ERR)
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
     1      call errquit('createfile: exceeded max_files',nfiles,
     2      UNKNOWN_ERR)
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
          call errquit('createfile: file problem',0,UNKNOWN_ERR)
  110     continue
          if (eaf_open(filename, EAF_RW, d_temp).ne.0)
     1      call errquit('createfile: eaf problem',0,DISK_ERR)
          if (eaf_close(d_temp).ne.0)
     1      call errquit('createfile: eaf problem',1,DISK_ERR)
          if (nodezero.and.util_print('files',print_debug)) then
            do i = 1,nfiles
              write(LuOut,9010) filehandles(i),filenames(i)(1:40)
            enddo
          endif
        else
          if (nfiles.eq.max_files)
     1      call errquit('createfile: exceeded max_files',nfiles,
     2      UNKNOWN_ERR)
          if (eaf_open(filename, EAF_RW, d_sf).ne.0)
     1      call errquit('createfile: eaf problem',0,DISK_ERR)
          nfiles = nfiles + 1
          filehandles(nfiles) = d_sf
          filenames(nfiles) = filename
        endif
        call tce_zero(d_sf,size)
      else if (ioalg.eq.2) then
        if (nodezero.and.util_print('files',print_debug))
     1    write(LuOut,9000) size,filename(1:40)
        if (.not.ga_create(mt_dbl,size,1,filename,
     1    -1,1,d_sf)) then
          write(LuOut,*) ' available GA memory ',
     1    ga_memory_avail(),' bytes'
          call errquit ('createfile: failed ga_create size/nproc bytes'
     S      ,   (size*ma_sizeof(mt_dbl,1,mt_byte))/ga_nnodes(),
     1      GA_ERR)
        endif
        call tce_zero(d_sf,size)
      else if (ioalg.eq.3) then
        if (parallel) call ga_sync()
        if (sf_create(filename,dble(bytes)*dble(size),
     1      dble(bytes)*dble(size),1.0d4,d_sf).ne.0)
     2    call errquit('createfile: sf problem',0,DISK_ERR)
        if (parallel) then
          if (sf_close(d_sf).ne.0)
     1      call errquit('createfile: sf problem',1,DISK_ERR)
          call ga_sync()
        endif
        if (nodezero.and.util_print('files',print_debug))
     1    write(LuOut,9010) d_sf,filename(1:40)
        call tce_zero(d_sf,size)
      else if (ioalg.eq.4) then
        if (nfiles.eq.max_files)
     1    call errquit('createfile: exceeded max_files',nfiles,
     2    UNKNOWN_ERR)
        if (eaf_open(filename, EAF_RW, d_sf).ne.0)
     1    call errquit('createfile: eaf problem',0,DISK_ERR)
        nfiles = nfiles + 1
        filehandles(nfiles) = d_sf
        filenames(nfiles) = filename
        if (util_print('files',print_debug))
     1    write(LuOut,9010) d_sf,filename(1:40)
        call tce_zero(d_sf,size)
      else if (ioalg.eq.5) then
        if (nfiles.eq.max_files)
     1    call errquit('createfile: exceeded max_files',nfiles,
     2    UNKNOWN_ERR)
        if (.not.ga_create(mt_dbl,size,1,filename,
     1    -1,1,d_temp))
     2    call errquit ('createfile: ga problem',1,GA_ERR)
        call ga_zero(d_temp)
        nfiles = nfiles + 1
        filehandles(nfiles) = d_temp
        filenames(nfiles) = filename
        filesizes(nfiles) = size
        fileisga(nfiles) = .true.
        d_sf = nfiles
      else if (ioalg.eq.6) then
        if (nfiles.eq.max_files)
     1    call errquit('createfile: exceeded max_files',nfiles,
     2    UNKNOWN_ERR)
        if (.not.ga_create(mt_dbl,size,1,filename,
     1    -1,1,d_temp))
     2    call errquit ('createfile: ga problem',1,GA_ERR)
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
















      subroutine createfile_sf(filename,d_sf,size)
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
      character*(*) filename
      integer d_sf
      integer d_temp
      integer size
      integer i
      logical nodezero
      logical parallel
      logical used
c
      nodezero = (ga_nodeid().eq.0)
      parallel = (ga_nnodes().gt.1)
      cpusecs(1) = cpusecs(1) - util_cpusec()
      cpusecs(51) = cpusecs(51) - util_wallsec()
      if (size.eq.0) then
        if (nodezero.and.util_print('files',print_high)) then
          write(LuOut,*) '*************** Warning ***************'
          write(LuOut,*) ' Create file size is zero. Calculation'
          write(LuOut,*) ' will continue by increasing the size.'
          write(LuOut,*) ' Use of a larger basis set is advised.'
        endif
        size = 1
      endif
c
c      else if (ioalg.eq.3) then
c
        if (parallel) call ga_sync()
        if (sf_create(filename,dble(bytes)*dble(size),
     1      dble(bytes)*dble(size),1.0d4,d_sf).ne.0)
     2    call errquit('createfile_sf: sf problem',0,DISK_ERR)
        if (parallel) then
          if (sf_close(d_sf).ne.0)
     1      call errquit('createfile_sf: sf problem',1,DISK_ERR)
          call ga_sync()
        endif
        if (nodezero.and.util_print('files',print_debug))
     1    write(LuOut,9010) d_sf,filename(1:40)
c
c        call tce_zero(d_sf,size)
c
        call tce_zero_sf(d_sf,size)
c
c      endif
c
      cpusecs(1) = cpusecs(1) + util_cpusec()
      cpusecs(51) = cpusecs(51) + util_wallsec()
 9000 format(1x,'create a file: size = ',i10,' file name = ',a)
 9010 format(1x,'file handle = ',i9,' file name = ',a)
 9020 format(1x,'node',i3,' file handle = ',i9,' file name = ',a)
 9030 format(1x,'create a file: node',i3,' size = ',i10,
     1          ' file name = ',a)
      return
      end












      subroutine createfile_sf_suffix(filename,d_sf,size,suffix)
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
      character*(*) filename
      integer d_sf
      integer d_temp
      integer size
      integer i
      integer suffix
      logical nodezero
      logical parallel
      logical used
c
      nodezero = (ga_nodeid().eq.0)
      parallel = (ga_nnodes().gt.1)
      cpusecs(1) = cpusecs(1) - util_cpusec()
      cpusecs(51) = cpusecs(51) - util_wallsec()
      if (size.eq.0) then
        if (nodezero.and.util_print('files',print_high)) then
          write(LuOut,*) '*************** Warning ***************'
          write(LuOut,*) ' Create file size is zero. Calculation'
          write(LuOut,*) ' will continue by increasing the size.'
          write(LuOut,*) ' Use of a larger basis set is advised.'
        endif
        size = 1
      endif
c
c      else if (ioalg.eq.3) then
c
        if (parallel) call ga_sync()
        if (sf_create_suffix(filename,dble(bytes)*dble(size),
     1      dble(bytes)*dble(size),1.0d4,d_sf,suffix).ne.0)
     2    call errquit('createfile_sf: sf problem',0,DISK_ERR)
        if (parallel) then
          if (sf_close(d_sf).ne.0)
     1      call errquit('createfile_sf: sf problem',1,DISK_ERR)
          call ga_sync()
        endif
        if (nodezero.and.util_print('files',print_debug))
     1    write(LuOut,9010) d_sf,filename(1:40)
c
c        call tce_zero(d_sf,size)
c
        call tce_zero_sf(d_sf,size)
c
c      endif
c
      cpusecs(1) = cpusecs(1) + util_cpusec()
      cpusecs(51) = cpusecs(51) + util_wallsec()
 9000 format(1x,'create a file: size = ',i10,' file name = ',a)
 9010 format(1x,'file handle = ',i9,' file name = ',a)
 9020 format(1x,'node',i3,' file handle = ',i9,' file name = ',a)
 9030 format(1x,'create a file: node',i3,' size = ',i10,
     1          ' file name = ',a)
      return
      end








      subroutine createfile_mirror(filename,d_sf,size)
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
      integer size
      integer proc_group
      logical nodezero
c      logical parallel
      logical status
      character*(*) filename
c
      nodezero = (ga_nodeid().eq.0)
c      parallel = (ga_nnodes().gt.1)
c
      cpusecs(1) = cpusecs(1) - util_cpusec()
      cpusecs(51) = cpusecs(51) - util_wallsec()
c
      if (size.eq.0) then
        if (nodezero.and.util_print('files',print_high)) then
          write(LuOut,*) '*************** Warning ***************'
          write(LuOut,*) ' Create file size is zero. Calculation'
          write(LuOut,*) ' will continue by increasing the size.'
          write(LuOut,*) ' Use of a larger basis set is advised.'
        endif
        size = 1
      endif
c
      if (nodezero.and.util_print('files',print_debug)) then
        write(LuOut,9000) size,filename(1:40)
      endif
c
c      status = ga_create(mt_dbl,size,1,filename,-1,1,d_sf)
c
      proc_group = ga_pgroup_get_mirror()
      status = nga_create_config(MT_DBL,     ! double precision elements
     1                           1,          ! a one-dimensional array
     2                           size,       ! size
     3                           filename,   ! name
     4                           -1,         ! chunking
     5                           proc_group, ! processor group handle
     6                           d_sf)       ! handle to create GA
c
      if (.not.status) then
        write(LuOut,*) ' available GA space ',ga_memory_avail(),' bytes'
        call errquit ('createfile: failed ga_create size=',size,GA_ERR)
      endif
c
      call ga_zero(d_sf)
c
      cpusecs(1) = cpusecs(1) + util_cpusec()
      cpusecs(51) = cpusecs(51) + util_wallsec()
c
 9000 format(1x,'create a file: size = ',i10,' file name = ',a)
 9010 format(1x,'file handle = ',i9,' file name = ',a)
 9020 format(1x,'node',i3,' file handle = ',i9,' file name = ',a)
 9030 format(1x,'create a file: node',i3,' size = ',i10,
     1          ' file name = ',a)
      return
      end






