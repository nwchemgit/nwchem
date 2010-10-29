      subroutine pre_rotate(o,p,angle,xw,mwm,mwa,nwm,nwa,xwc,mwmc,nwmc,
     + xs,msa,nsa)
c
c $Id$
c
      implicit none
c
      integer mwm,mwa,nwm,nwa,mwmc,nwmc,msa,nsa
      real*8 xw(3,mwa,mwm),xwc(3,mwa,mwmc),xs(3,msa)
      real*8 angle,o(3),p(3),x(3),y(3)
c
      integer i,j
c
c     rotate
c
      do 1 i=1,nsa
      x(1)=xs(1,i)
      x(2)=xs(2,i)
      x(3)=xs(3,i)
      call rotate(o,p,angle,x,y)
      xs(1,i)=y(1)
      xs(2,i)=y(2)
      xs(3,i)=y(3)
    1 continue
      do 2 i=1,nwmc
      do 3 j=1,nwa
      x(1)=xwc(1,j,i)
      x(2)=xwc(2,j,i)
      x(3)=xwc(3,j,i)
      call rotate(o,p,angle,x,y)
      xwc(1,j,i)=y(1)
      xwc(2,j,i)=y(2)
      xwc(3,j,i)=y(3)
    3 continue
    2 continue
      do 4 i=1,nwm
      do 5 j=1,nwa
      x(1)=xw(1,j,i)
      x(2)=xw(2,j,i)
      x(3)=xw(3,j,i)
      call rotate(o,p,angle,x,y)
      xw(1,j,i)=y(1)
      xw(2,j,i)=y(2)
      xw(3,j,i)=y(3)
    5 continue
    4 continue
      return
      end
