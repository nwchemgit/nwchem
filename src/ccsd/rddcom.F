      subroutine rddcom()
      implicit integer (a-z)
#include "global.fh"
      common/debug/doa(15),dob(15),dog(15),doh(15),dojk(15),
     &             dos(15),dod(15)
c
c      open(unit=56,name='debug',status='old')
c      read(56,*)doa
c      read(56,*)dob
c      read(56,*)dog
c      read(56,*)doh
c      read(56,*)dojk
c      read(56,*)dos
c      read(56,*)dod
c      close(unit=56)
      do i=1,15
       doa(i)=2
       dob(i)=2
       dog(i)=1
       doh(i)=1
       dojk(i)=2
       dos(i)=1
       dod(i)=1
      enddo
c
      call ga_sync()
      if (ga_nodeid().eq.0)then
      write(*,57)' doa  ',doa
      write(*,57)' dob  ',dob
      write(*,57)' dog  ',dog
      write(*,57)' doh  ',doh
      write(*,57)' dojk ',dojk
      write(*,57)' dos  ',dos
      write(*,57)' dod  ',dod
   57 format(1x,a6,1x,15i2)
      endif
      call ga_sync()
c
      return
      end

