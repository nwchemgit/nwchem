      subroutine argos_cafe_pmeset(xw,yw,nwmloc,xs,ys,nsaloc,
     + theta,dtheta)
c
      implicit none
c
#include "argos_cafe_common.fh"
c
      real*8 xw(mwm,3,mwa),xs(msa,3)
      real*8 yw(mwm,3,mwa),ys(msa,3)
      real*8 theta(morder,3,*),dtheta(morder,3,*)
      integer nwmloc,nsaloc
c
      integer i,j,iwa,ix,iwm,isa,indx,ng(3)
      real*8 w,rinv
c
      if(nbxtyp.ne.1) then
      do 1 i=1,3
      do 2 j=1,3
      recip(i,j)=zero
    2 continue
      recip(i,i)=one/box(i)
    1 continue
      else
      do 3 i=1,3
      do 4 j=1,3
      recip(i,j)=vlati(i,j)
    4 continue
    3 continue
      endif
c
c     calculate scaled fractional coordinates
c
      ng(1)=ngx
      ng(2)=ngy
      ng(3)=ngz
      do 5 iwa=1,mwa
      do 6 ix=1,3
      do 7 iwm=1,nwmloc
      w=xw(iwm,1,iwa)*recip(1,ix)+xw(iwm,2,iwa)*recip(2,ix)+
     + xw(iwm,3,iwa)*recip(3,ix)
      yw(iwm,ix,iwa)=dble(ng(ix))*(w-anint(w)+half)
    7 continue
    6 continue
    5 continue
c
      do 8 ix=1,3
      do 9 isa=1,nsaloc
      w=xs(isa,1)*recip(1,ix)+xs(isa,2)*recip(2,ix)+
     + xs(isa,3)*recip(3,ix)
      ys(isa,ix)=dble(ng(ix))*(w-anint(w)+half)
    9 continue
    8 continue
c
c     spline coefficients
c
      indx=0
      do 10 iwa=1,mwa
      do 11 iwm=1,nwmloc
      indx=indx+1
      do 12 ix=1,3
      w=yw(iwm,ix,iwa)-int(yw(iwm,ix,iwa))
      theta(morder,ix,indx)=zero
      theta(2,ix,indx)=w
      theta(1,ix,indx)=one-w
      do 13 i=3,morder-1
      rinv=one/dble(i-1)
      theta(i,ix,indx)=rinv*w*theta(i-1,ix,indx)
      do 14 j=1,i-2
      theta(i-j,ix,indx)=rinv*((w+dble(j))*theta(i-j-1,ix,indx)+
     + (dble(i-j)-w)*theta(i-j,ix,indx))
   14 continue
      theta(1,ix,indx)=rinv*(one-w)*theta(1,ix,indx)
   13 continue
      dtheta(1,ix,indx)=-theta(1,ix,indx)
      do 15 j=2,morder
      dtheta(j,ix,indx)=theta(j-1,ix,indx)-theta(j,ix,indx)
   15 continue
      rinv=one/dble(morder-1)
      theta(morder,ix,indx)=rinv*w*theta(morder-1,ix,indx)
      do 16 j=1,morder-2
      theta(morder-j,ix,indx)=rinv*((w+dble(j))*
     + theta(morder-j-1,ix,indx)+
     + (dble(morder-j)-w)*theta(morder-j,ix,indx))
   16 continue
      theta(1,ix,indx)=rinv*(one-w)*theta(1,ix,indx)
   12 continue
   11 continue
   10 continue
c
      do 17 isa=1,nsaloc
      indx=indx+1
      do 18 ix=1,3
      w=ys(isa,ix)-int(ys(isa,ix))
      theta(morder,ix,indx)=zero
      theta(2,ix,indx)=w
      theta(1,ix,indx)=one-w
      do 19 i=3,morder-1
      rinv=one/dble(i-1)
      theta(i,ix,indx)=rinv*w*theta(i-1,ix,indx)
      do 20 j=1,i-2
      theta(i-j,ix,indx)=rinv*((w+dble(j))*theta(i-j-1,ix,indx)+
     + (dble(i-j)-w)*theta(i-j,ix,indx))
   20 continue
      theta(1,ix,indx)=rinv*(one-w)*theta(1,ix,indx)
   19 continue
      dtheta(1,ix,indx)=-theta(1,ix,indx)
      do 21 j=2,morder
      dtheta(j,ix,indx)=theta(j-1,ix,indx)-theta(j,ix,indx)
   21 continue
      rinv=one/dble(morder-1)
      theta(morder,ix,indx)=rinv*w*theta(morder-1,ix,indx)
      do 22 j=1,morder-2
      theta(morder-j,ix,indx)=rinv*((w+dble(j))*
     + theta(morder-j-1,ix,indx)+
     + (dble(morder-j)-w)*theta(morder-j,ix,indx))
   22 continue
      theta(1,ix,indx)=rinv*(one-w)*theta(1,ix,indx)
   18 continue
   17 continue
c
      return
      end
c $Id$
