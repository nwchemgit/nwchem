      Subroutine dft_invio(name,iga_inv,nbf,ioop,lmiss)

C     $Id: dft_invio.F,v 1.7 1996-10-18 21:50:12 d3h449 Exp $
      implicit none
      
      character*(*) name        ! [input] name of file
      character*(*) ioop        ! [input] operation READ/WRITE/KILL
      integer lmiss             ! [output] 0=success, 1=failure
      integer iunit,j,nbf,iga_inv,istart,i
      integer iinv,linv
#include "mafdecls.fh"
#include "stdio.fh"
#include "util.fh"
      character*(nw_max_path_len) fname
c     
c     Hacked by RJH to close files after all operations and to 
c     use only one unit number ... old values left in comments
c     
      iunit = 77
      if (name.eq.'XC') then
         call util_file_name('xc', .true., .false., fname)
*         iunit=78
      elseif (name.eq.'CD') then
         call util_file_name('cdfit', .true., .false., fname)
*         iunit=76
      elseif (name.eq.'CDI') then
         call util_file_name('cdinv', .true., .false., fname)
*         iunit=77
      elseif (name.eq.'FK') then
         call util_file_name('fock', .true., .false., fname)
*         iunit=20
      endif
c     
      if(ioop.eq.'KILL') then
         call util_file_unlink(fname)
         lmiss = 0
         return
      endif
c     
      if(.not.ma_push_get(MT_Dbl,nbf,'inv i/o',linv,iinv))
     &     call errquit('dft_invio: cannot allocate inv i/o',0)
c     
      open(unit=iunit,file=fname,access='sequential',
     &     form='unformatted',status='unknown')
c     
      if(ioop.eq.'READ') then
         istart=iinv
         do j=1,nbf
            read(iunit,err=27,end=27)(Dbl_MB(iinv+i-1),i=1,
     &           nbf)
            call put_col(iga_inv,nbf,j,Dbl_MB(istart))
         enddo
         write(LuOut,*) ' READ file ', fname
         close(iunit)
         lmiss=0
         if(.not.ma_pop_stack(linv))
     &        call errquit('dft_invio: cannot pop stack',0)
         return
 27      lmiss=1
         if(.not.ma_pop_stack(linv))
     &        call errquit('dft_invio: cannot pop stack',0)
         close(iunit)
         return
      else
         istart=iinv
         lmiss=0
         do j=1,nbf
            call get_col(iga_inv,nbf,j,Dbl_MB(istart))
            write(iunit)(Dbl_MB(iinv+i-1),i=1,nbf)
         enddo
         if(.not.ma_pop_stack(linv))
     &        call errquit('dft_invio: cannot pop stack',0)
         close(iunit)
         return
      endif
      END
