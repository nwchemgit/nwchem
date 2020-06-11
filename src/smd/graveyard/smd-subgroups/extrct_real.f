c
c $Id$
c

      SUBROUTINE ex_1real(n,ist,irecord,buf,ierr)

      implicit none

      integer i,j,n,ifact2,ifact,ist,nstart
      integer ierr,limit,iexft,irec_len

      real*8 expo,buf,fact

      character*1 irecord,cchar,itemp

      parameter(irec_len=100)

      dimension irecord(irec_len),cchar(19)

      data cchar/'0','1','2','3','4','5','6','7','8','9','+','&','/',
     $           '-','.','D','E','d','e'/

      buf=0.0d+0
      expo=0.0d+0
      iexft=0
      limit=19
      ifact2=0
      fact=1.0d0
      nstart=ist+n-1
      do i=1,n
       do j=1,limit
        if(cchar(j).eq.irecord(nstart))goto 180
       enddo
 170   ierr=2
       return
 180   if(j.lt.11)goto 200
       if(j.le.14)goto 190
       if(j.gt.15)then
        expo=buf
        iexft=i
        ifact2=i
        buf=0.0d+0
        fact=1.0d0
        goto 210
       endif
       ifact2=i-1
       limit=14
       goto 210
 190   continue
       if(j.eq.14)buf=-buf
       goto 210
 200   buf=buf+(dble(j-1)*fact)
       fact=fact*10.0d+0
 210   nstart=nstart-1
      enddo
 220  buf=(0.1d0**(ifact2-iexft))*buf
      buf=buf*10**expo
 230  continue

      return

      END
