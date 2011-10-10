      subroutine argos_addref(ladd,xwm,xwcr,xsm,xscr,dsr)
c
      implicit none
c
#include "argos.fh"
c
      logical ladd
      real*8 xwm(mwm,3),xwcr(mwm,3),xsm(msm,3),xscr(msm,3),dsr(msm)
c
      integer i
c
c      return
      if(ladd) then
      dwr=zero
      do 1 i=1,nwmloc
      xwcr(i,1)=xwcr(i,1)+xwm(i,1)
      xwcr(i,2)=xwcr(i,2)+xwm(i,2)
      xwcr(i,3)=xwcr(i,3)+xwm(i,3)
      dwr=dwr+xwcr(i,1)*xwcr(i,1)+xwcr(i,2)*xwcr(i,2)+
     + xwcr(i,3)*xwcr(i,3)
    1 continue
      do 2 i=1,nsm
      xscr(i,1)=xscr(i,1)+xsm(i,1)
      xscr(i,2)=xscr(i,2)+xsm(i,2)
      xscr(i,3)=xscr(i,3)+xsm(i,3)
      dsr(i)=xscr(i,1)*xscr(i,1)+xscr(i,2)*xscr(i,2)+xscr(i,3)*xscr(i,3)
    2 continue
      else
      do 3 i=1,nwmloc
      xwcr(i,1)=xwcr(i,1)-xwm(i,1)
      xwcr(i,2)=xwcr(i,2)-xwm(i,2)
      xwcr(i,3)=xwcr(i,3)-xwm(i,3)
    3 continue
      do 4 i=1,nsm
      xscr(i,1)=xscr(i,1)-xsm(i,1)
      xscr(i,2)=xscr(i,2)-xsm(i,2)
      xscr(i,3)=xscr(i,3)-xsm(i,3)
    4 continue
      endif
c
      return
      end
c $Id$
