      Subroutine dft_invio(name,iga_inv,nbf,ioop,lmiss)

C     $Id: dft_invio.F,v 1.3 1995-10-17 08:55:12 d3g681 Exp $
      Implicit real*8 (a-h,o-z)
      Implicit integer (i-n)
      character*10 fname
      character*(*) name,ioop
      logical LResult, doesit, isit
#include "mafdecls.fh"
#include "stdio.fh"
      
      if (name.eq.'XC') then
        fname='xc_inv.dat'
        iunit=78
      elseif (name.eq.'CD') then
        fname='cd.dat'
        iunit=76
      elseif (name.eq.'CDI') then
        fname='cd_inv.dat'
        iunit=77
      elseif (name.eq.'FK') then
        fname='fock.dat'
        iunit=20
      endif
      
      inquire(file=fname,exist=doesit,opened=isit)
      if(ioop.eq.'KILL'.and.doesit) then
        if(.not.isit)
     &       open(unit=iunit,file=fname,access='sequential',
     &       form='unformatted',status='old')
        close(iunit,status='delete')
        return
      endif
      LResult = MA_Push_Get(MT_Dbl,nbf,'inv i/o',linv,iinv)
      
      if(.not.isit)then
        if(.not.doesit)then
          open(unit=iunit,file=fname,access='sequential',
     &         form='unformatted',status='new')
        else
          open(unit=iunit,file=fname,access='sequential',
     &         form='unformatted',status='old')
        endif
      endif

      if(ioop.eq.'READ') then
        istart=iinv
        do j=1,nbf
          read(iunit,err=27,end=27)(Dbl_MB(iinv+i-1),i=1,
     &         nbf)
          call put_col(iga_inv,nbf,j,Dbl_MB(istart))
        enddo
        write(LuOut,*) ' READ file ',fname
        rewind iunit
        lmiss=0
        LResult = MA_Pop_Stack(linv)
        close(iunit)
        return
   27   lmiss=1
        LResult = MA_Pop_Stack(linv)
        rewind iunit
C       close(iunit)
        return
      else
        istart=iinv
        lmiss=0
        do j=1,nbf
          call get_col(iga_inv,nbf,j,Dbl_MB(istart))
          write(iunit)(Dbl_MB(iinv+i-1),i=1,nbf)
        enddo
        LResult = MA_Pop_Stack(linv)
        rewind iunit
C       close(iunit)
        return
      endif
      END
