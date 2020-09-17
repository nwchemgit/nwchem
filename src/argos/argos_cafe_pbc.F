      subroutine argos_cafe_pbc(iopt,xd,ndd,xc,ndc,noff,n0,n1)
c
c $Id$
c
      implicit none
c
#include "argos_cafe_common.fh"
c
      integer iopt,ndd,ndc,noff,n0,n1
      real*8 xd(ndd,3),xc(ndc,3),xt(3),xts
      integer i,ix,jx,ndim
      real*8 rin(3),pbcor,boxi(3)
c
      if(npbtyp.eq.0) return
c
      boxi(1)=one/box(1)
      boxi(2)=one/box(2)
      boxi(3)=one/box(3)
c
c     correction vector returned in xc
c     --------------------------------
c
      if(iopt.eq.0) then
c
      do 1 i=n0,n1
      xc(i,1)=zero
      xc(i,2)=zero
      xc(i,3)=zero
    1 continue
c
      ndim=3
      if(npbtyp.eq.2) ndim=2
c
c     rectangular box
c
      if(nbxtyp.eq.0.or.nbxtyp.eq.2) then
      do 2 ix=1,ndim
      do 3 i=n0,n1
      if(abs(xd(noff+i,ix)).gt.boxh(ix)) then
      xc(i,ix)=anint(xd(noff+i,ix)*boxi(ix))*box(ix)
      endif
    3 continue
    2 continue
      endif
c
      if(nbxtyp.eq.1) then
      do 4 i=n0,n1
      do 5 ix=1,3
      xts=zero
      do 6 jx=1,3
      xts=xts+vlati(ix,jx)*xd(noff+i,jx)
    6 continue
      xt(ix)=int(two*xts)
    5 continue
      do 7 ix=1,3
      xc(i,ix)=zero
      do 8 jx=1,3
      xc(i,ix)=xc(i,ix)+vlat(ix,jx)*xt(jx)
    8 continue
    7 continue
    4 continue
      endif
c
      if(nbxtyp.eq.2) then
      do 9 i=n0,n1
      pbcor=boxh(1)*aint(four*(abs((xd(noff+i,1)-xc(i,1))/box(1))
     + +abs((xd(noff+i,2)-xc(i,2))/box(2))
     + +abs((xd(noff+i,3)-xc(i,3))/box(3)))/three)
      do 10 ix=1,3
      xc(i,ix)=xc(i,ix)+sign(pbcor,(xd(noff+i,ix)-xc(i,ix)))
   10 continue
    9 continue
      endif
c
      endif
c
c     correction applied to xd
c     ------------------------
c
      if(iopt.eq.1) then
c
      ndim=3
      if(npbtyp.eq.2) ndim=2
c
c     rectangular box
c
      if(nbxtyp.eq.0.or.nbxtyp.eq.2) then
      do 12 ix=1,ndim
      do 13 i=n0,n1
      if(abs(xd(noff+i,ix)).gt.boxh(ix))
     + xd(noff+i,ix)=xd(noff+i,ix)-anint(xd(noff+i,ix)*boxi(ix))*box(ix)
   13 continue
   12 continue
      endif
c
      if(nbxtyp.eq.1) then
      do 14 i=n0,n1
      do 15 ix=1,3
      rin(ix)=xd(noff+i,ix)
      xts=zero
      do 16 jx=1,3
      xts=xts+vlati(ix,jx)*xd(noff+i,jx)
   16 continue
      xt(ix)=int(two*xts)
   15 continue
      do 17 ix=1,3
      do 18 jx=1,3
      xd(noff+i,ix)=xd(noff+i,ix)-vlat(ix,jx)*xt(jx)
   18 continue
   17 continue
   14 continue
      endif
c
      if(nbxtyp.eq.2) then
      do 19 i=n0,n1
      pbcor=boxh(1)*aint(four*(abs(xd(noff+i,1)/box(1))
     + +abs(xd(noff+i,2)/box(2))+abs(xd(noff+i,3)/box(3)))/three)
      do 20 ix=1,3
      xd(noff+i,ix)=xd(noff+i,ix)-sign(pbcor,xd(noff+i,ix))
   20 continue
   19 continue
      endif
c
      endif
c
      return
      end
