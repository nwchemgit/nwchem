      subroutine dia_xaver(xtrj,ntrj,nfr,xa,nxa)
c
c      x,xa,nx,mx,nave)
c
c $Id: dia_xaver.F,v 1.2 2007-01-03 05:43:27 d3j191 Exp $
c
      implicit none
c
#include "msgids.fh"
#include "global.fh"
c
      real*8 zero
      parameter(zero=0.0d0)
c
      integer ntrj,nfr,nxa
      real*8 xtrj(ntrj,3,nfr),xa(nxa,3)
c
c      integer nx,mx,nave
c      real*8 x(mx,3),xa(mx,3)
c
      integer i,j,k
c
      do 1 j=1,3
      do 2 i=1,nxa
      xa(i,j)=zero
    2 continue
      do 3 k=1,nfr
      do 4 i=1,ntrj
      xa(i,j)=xa(i,j)+xtrj(i,j,k)
    4 continue
    3 continue
    1 continue
c
c     global sum
c
      call ga_dgop(mag_d08,xa,3*nxa,'+')
c
      write(*,'(a)') 'SUM'
      write(*,'(i5,3f12.6)') (i,(xa(i,j),j=1,3),i=1,10)
c
c      if(nave.eq.0) then
c      do 1 j=1,3
c      do 2 i=1,nx
c      xa(i,j)=x(i,j)
c    2 continue
c    1 continue
c      else
c      do 3 j=1,3
c      do 4 i=1,nx
c      xa(i,j)=xa(i,j)+x(i,j)
c    4 continue
c    3 continue
c      endif
c
c      nave=nave+1
c
      return
      end
      subroutine dia_xmean(xa,nx,mx,nave)
c
      implicit none
c
#include "msgids.fh"
#include "global.fh"
c
      integer nx,mx,nave
      real*8 xa(mx,3)
      integer i,j
      real*8 fact
c
      fact=1.0d0/dble(nave)
c
      do 1 j=1,3
      do 2 i=1,nx
      xa(i,j)=xa(i,j)*fact
    2 continue
    1 continue
c
      return
      end
