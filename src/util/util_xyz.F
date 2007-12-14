      function util_xyz_nframes(in_xyz,nf)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "stdio.fh"
#include "inp.fh"
      character*(*) in_xyz
      integer nf
      logical util_xyz_nframes
c     local variables
      character*32 pname
      integer i,ns
      character*16 buf
      character*255 filename
      character*255 xyzfile
      character*255 message

      integer fn_xyz

      logical util_get_io_unit
      external util_get_io_unit

      pname = "util_xyz_nframes"
      util_xyz_nframes = .false.
c
      xyzfile = in_xyz(1:inp_strlen(in_xyz))
      call util_file_name_resolve(xyzfile, .false.)
c
      filename = in_xyz(1:inp_strlen(xyzfile))
c
c     prepare files for reading/writing
c     ---------------------------------
      if(.not.util_get_io_unit(fn_xyz)) 
     >       call errquit("cannot get file number",0,0)
      filename = xyzfile
      open(fn_xyz,file=filename,form='formatted',status='old',
     $          err=133)
c
      nf = 0
c
c     get number of atoms
c     ------------------
  9   continue
      message = " number of atoms "
      read(fn_xyz,*,err=10,end=10) ns
c
c     seek frame
c     ----------
      message = " seek frame"
      do i=1,ns+1
        read(fn_xyz,*,err=10,end=10) buf 
c        for strict/hard processing use
c        read(fn_xyz,*,err=134,end=135) buf 
      end do
      nf = nf + 1
      goto 9
 10   continue 
      close(fn_xyz)
      util_xyz_nframes = .true.
      return
 133  call errquit(pname//'error opening/closing '//filename,0, 0)
 134  call errquit(pname//'error reading xyz file'//message,0, 0)
 135  call errquit(pname//'error end of file at'//message,0, 0)

      end

      function util_xyz_seek(fn_xyz,nf)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "stdio.fh"
#include "inp.fh"
      integer fn_xyz
      integer nf
      logical util_xyz_seek
c     local variables
      integer sf
      character*32 pname
      integer i,ns
      character*16 buf
      character*255 message

      pname = "util_xyz_seek"
      util_xyz_seek = .false.
c
c     prepare files for reading/writing
c     ---------------------------------
      rewind(fn_xyz)
      sf = 0
c
c     get number of atoms
c     ------------------
  9   continue
      message = " number of atoms "
      read(fn_xyz,*,err=10,end=10) ns
c
c     seek frame
c     ----------
      message = " seek frame"
      do i=1,ns+1
        read(fn_xyz,*,err=10,end=10) buf 
      end do
      sf = sf + 1
      if(sf.eq.(nf-1)) then
        util_xyz_seek = .true.
        return
      end if 
      goto 9
 10   continue 
      return

      end
