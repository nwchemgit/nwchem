      logical function pre_repeat(card,drep,
     + csa,isat,isgm,xs,vs,msa,nsa,xw,vw,mwm,mwa,nwm,nwa,
     + xwc,vwc,mwmc,nwmc)
c
c $Id: pre_repeat.F,v 1.1 1999-11-03 15:43:51 d3j191 Exp $
c
      implicit none
c
#include "util.fh"
c
      integer msa,nsa,mwm,nwm,mwa,nwa,mat,mwmc,nwmc
      integer isat(msa),isgm(msa)
      character*16 csa(msa)
      real*8 xs(3,msa),vs(3,msa),xw(3,mwa,mwm),vw(3,mwa,mwm)
      real*8 xwc(3,mwa,mwmc),vwc(3,mwa,mwmc),drep(3)
      character*80 card
c
      integer i,j,k,imin(3),jmin(3)
      real*8 xmax,xmin,r2min,dx(3),ds(3),r2min(3),r2(3),dd(3,3)
      real*8 touch
      logical done
c
      read(card(23:34),'(f12.6)') touch
c
      do 1 k=1,3
      xmax=xs(k,1)
      xmin=xs(k,1)
      do 2 i=2,nsa
      if(xs(k,i).gt.xmax) xmax=xs(k,i)
      if(xs(k,i).lt.xmin) xmin=xs(k,i)
    2 continue
      dx(k)=xmax-xmin+touch
    1 continue
c
   10 continue
c 
      do 3 k=1,3
      imin(k)=1
      jmin(k)=1
      r2min(k)=1.0e24
    3 continue
c
      do 4 i=1,nsa
      do 5 j=1,nsa
      do 6 k=1,3
      ds(k)=xs(k,i)-xs(k,j)
    6 continue
      r2(1)=(ds(1)-dx(1))**2+ds(2)**2+ds(3)**2
      r2(2)=ds(1)**2+(ds(2)-dx(2))**2+ds(3)**2
      r2(3)=ds(1)**2+ds(2)**2+(ds(3)-dx(3))**2
      do 7 k=1,3
      if(r2(k).lt.r2min(k)) then
      imin(k)=i
      jmin(k)=j
      r2min(k)=r2(k)
      endif
    7 continue
    5 continue
    4 continue
c
      write(*,'(a,3f12.6)') 'x ',(dx(k),k=1,3)
      write(*,'(a,3f12.6)') 'r ',(sqrt(r2min(k)),k=1,3)
c
      done=.true.
c
      do 9 k=1,3
      if(sqrt(r2min(k)).gt.2.0d0*touch) then
      done=.false.
      dx(k)=dx(k)-touch
      endif
    9 continue
c
      if(.not.done) goto 10
c
      do 8 k=1,3
      drep(k)=dx(k)+touch-
     + abs(xs(k,imin(k))-xs(k,jmin(k))-dx(k))
    8 continue
      write(*,'(a,3f12.6)') 'd ',(drep(k),k=1,3)
c
      pre_repeat=.true.
      return
c
 9999 continue
      pre_repeat=.false.
      return
      end
