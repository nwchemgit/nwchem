      subroutine dia_xaver(xtrj,ntrj,nfr,xa,nxa)
c
c      x,xa,nx,mx,nave)
c
c $Id: dia_xaver.F,v 1.4 2007-01-12 23:43:02 d3j191 Exp $
c
      implicit none
c
#include "msgids.fh"
#include "global.fh"
#include "dia_common.fh"
c
      real*8 zero
      parameter(zero=0.0d0)
c
      integer ntrj,nfr,nxa
      real*8 xtrj(ntrj,3,nfr),xa(nxa,3)
c
      integer i,j,k
c
      if(nxa.ne.ntrj) call md_abort('xaver dimension problem',0)
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
      do 5 j=1,3
      do 6 i=1,ntrj
      xa(i,j)=xa(i,j)/dble(nfrtot)
    6 continue
    5 continue
c
ctest
c      if(me.eq.0) then
c      do 7 i=1,nfr
c      write(*,'(i5,3f12.6)') i,(xtrj(5,j,i),j=1,3)
c    7 continue
c      endif
ctest
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
