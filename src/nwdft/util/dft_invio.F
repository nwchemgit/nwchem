      Subroutine dft_invio(name,iga_inv,nbf,ioop,lmiss)

C     $Id$
      implicit none
#include "errquit.fh"
      
      character*(*) name        ! [input] name of file
      character*(*) ioop        ! [input] operation READ/WRITE/KILL
      integer lmiss             ! [output] 0=success, 1=failure
      integer j,nbf,iga_inv,i
      integer iinv,linv,me,ierr,fd
#include "mafdecls.fh"
#include "stdio.fh"
#include "util.fh"
#include "eaf.fh"
#include "global.fh"
      character*(nw_max_path_len) fname
      logical oprint_ioinfo
      double precision offs
c     
c     Hacked by RJH to close files after all operations and to 
c     use only one unit number ... old values left in comments
c     
      me=ga_nodeid()
      oprint_ioinfo = util_print('io info',print_debug)
      if (name.eq.'XC') then
         call util_file_name('xc', .true., .false., fname)
      elseif (name.eq.'CD') then
         call util_file_name('2ceri', .true., .false., fname)
      elseif (name.eq.'CDI') then
         call util_file_name('2ceri_inv', .true., .false., fname)
      elseif (name.eq.'FK') then
         call util_file_name('fock', .true., .false., fname)
      endif
c     
      if(ioop.eq.'KILL') then
         call util_file_unlink(fname)
         ierr=eaf_delete(fname)
         if(ierr.ne.0) call errquit('dftinvio: eaf_delete',ierr,
     &       DISK_ERR)
         if(me.eq.0.and.oprint_ioinfo)
     &   write(LuOut,*) ' KILL file ', fname, ' lmiss = ',lmiss
         lmiss = 0
         return
      endif
c     
      if(.not.ma_push_get(MT_Dbl,nbf,'inv i/o',linv,iinv))
     &     call errquit('dft_invio: cannot allocate inv i/o',0, MA_ERR)
c     
      ierr=eaf_open(fname, eaf_rw, fd)
      if(ierr.ne.0) call errquit('dftinvio: eaf_open',ierr, DISK_ERR)
c     
      offs=0d0
      if(ioop.eq.'READ') then
         do j=1,nbf
            ierr = eaf_read(fd, offs, dbl_mb(iinv),8*nbf)
      if(ierr.ne.0) call errquit('dftinvio: eaf_read',ierr,
     &       DISK_ERR)
            offs=offs+8*nbf
            call put_col(iga_inv,nbf,j,Dbl_MB(iinv))
         enddo
         ierr=eaf_close(fd)
      if(ierr.ne.0) call errquit('dftinvio: eaf_close1',ierr,
     &       DISK_ERR)
         lmiss=0
         if(.not.ma_pop_stack(linv))
     &        call errquit('dft_invio: cannot pop stack',0, MA_ERR)
         if(me.eq.0.and.oprint_ioinfo)
     &   write(LuOut,*) ' READ file ', fname, ' lmiss = ',lmiss
         return
 27      lmiss=1
         if(.not.ma_pop_stack(linv))
     &        call errquit('dft_invio: cannot pop stack',0, MA_ERR)
         ierr=eaf_close(fd)
         if(ierr.ne.0) call errquit('dftinvio: eaf_close2',ierr,
     &       DISK_ERR)
         if(me.eq.0.and.oprint_ioinfo)
     &   write(LuOut,*) ' READ file ', fname, ' lmiss = ',lmiss
         return
      else
         lmiss=0
         offs=0d0
         do j=1,nbf
            call get_col(iga_inv,nbf,j,Dbl_MB(iinv))
            ierr = eaf_write(fd, offs, dbl_mb(iinv),8*nbf)
            if(ierr.ne.0) call errquit('dftinvio: eaf_write',ierr,
     &       DISK_ERR)
            offs=offs+8*nbf
         enddo
         if(.not.ma_pop_stack(linv))
     &        call errquit('dft_invio: cannot pop stack',0, MA_ERR)
         ierr=eaf_close(fd)
         if(ierr.ne.0) call errquit('dftinvio: eaf_close3',ierr,
     &       DISK_ERR)
         if(me.eq.0.and.oprint_ioinfo)
     &   write(LuOut,*) ' WRITE file ', fname, ' lmiss = ',lmiss
         return
      endif
      END
