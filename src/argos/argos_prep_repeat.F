      logical function argos_prep_repeat(card,lfnout,drep,box,
     + csa,isat,isgm,xs,vs,msa,nsa,xw,vw,mwm,mwa,nwm,nwa,
     + xwc,vwc,mwmc,nwmc)
c
c $Id$
c
      implicit none
c
#include "util.fh"
c
      integer msa,nsa,mwm,nwm,mwa,nwa,mwmc,nwmc,lfnout
      integer isat(msa),isgm(msa)
      character*16 csa(msa)
      real*8 xs(3,msa),vs(3,msa),xw(3,mwa,mwm),vw(3,mwa,mwm)
      real*8 xwc(3,mwa,mwmc),vwc(3,mwa,mwmc),drep(3),box(3)
      character*80 card
c
      integer i,j,k,imin(3),jmin(3)
      real*8 xmax,xmin,dx(3),ds(3),r2min(3),r2(3)
      real*8 touch,ddx(3)
      logical done
c
      read(card(23:34),'(f12.6)') touch
c
      if(touch.le.0.0d0.and.box(1).gt.0.0d0.and.
     + box(2).gt.0.0d0.and.box(3).gt.0.0d0) then
      drep(1)=box(1)
      drep(2)=box(2)
      drep(3)=box(3)
      argos_prep_repeat=.true.
      return
      endif
c
      do 1 k=1,3
      ddx(k)=touch
      xmax=xs(k,1)
      xmin=xs(k,1)
      do 2 i=2,nsa
      if(xs(k,i).gt.xmax) xmax=xs(k,i)
      if(xs(k,i).lt.xmin) xmin=xs(k,i)
    2 continue
      dx(k)=xmax-xmin+1.5d0*touch
    1 continue
c
      if(util_print('restart',print_high)) then
      write(lfnout,1000)
 1000 format(/,' Convergence',/,15x,'x',11x,'y',11x,'z',/)
      endif
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
      ds(k)=xs(k,j)-xs(k,i)
    6 continue
      r2(1)=(ds(1)+dx(1))**2+ds(2)**2+ds(3)**2
      r2(2)=ds(1)**2+(ds(2)+dx(2))**2+ds(3)**2
      r2(3)=ds(1)**2+ds(2)**2+(ds(3)+dx(3))**2
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
c      write(*,'(a,3f12.6)') 'x ',(dx(k),k=1,3)
c      write(*,'(a,3f12.6)') 'r ',(sqrt(r2min(k)),k=1,3)
c
      done=.true.
c
      do 9 k=1,3
      if(sqrt(r2min(k)).lt.touch) then
      dx(k)=dx(k)+ddx(k)
      ddx(k)=0.5d0*ddx(k)
      endif
      dx(k)=dx(k)-ddx(k)
      if(ddx(k).gt.0.01) done=.false.
    9 continue
c
      if(util_print('restart',print_high)) then
      write(lfnout,1001) ddx
 1001 format(10x,3f12.6)
      endif
c
      if(.not.done) goto 10
c
      do 8 k=1,3
      drep(k)=dx(k)
    8 continue
c
      if(util_print('restart',print_medium)) then
      write(lfnout,1002)
 1002 format(' Repeat separation determined')
      endif
c
      argos_prep_repeat=.true.
      return
c
 9999 continue
      argos_prep_repeat=.false.
      return
      end
