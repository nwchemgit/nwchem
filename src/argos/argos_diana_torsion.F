      real*8 function argos_diana_torsion(x,nx,w,nw,na,i,i2,j,j2,k,k2,l,
     + l2)
c
c $Id$
c
      implicit none
#include "global.fh"
c
      integer nx,nw,na,i,i2,j,j2,k,k2,l,l2
      real*8 x(nx,3),w(nw,na,3)
c
      real*8 xijx,xkjx,xklx,xikx,xjlx
      real*8 xijy,xkjy,xkly,xiky,xjly
      real*8 xijz,xkjz,xklz,xikz,xjlz
      real*8 xmx,xmy,xmz,xnx,xny,xnz,rm2i,rn2i,rmni,cphi,s,phi
c
      real*8 xa(3),xb(3),xc(3),xd(3)
c
      if(i2.le.0) then
      xa(1)=x(i,1)
      xa(2)=x(i,2)
      xa(3)=x(i,3)
      else
      xa(1)=w(i,i2,1)
      xa(2)=w(i,i2,2)
      xa(3)=w(i,i2,3)
      endif
      if(j2.le.0) then
      xb(1)=x(j,1)
      xb(2)=x(j,2)
      xb(3)=x(j,3)
      else
      xb(1)=w(j,j2,1)
      xb(2)=w(j,j2,2)
      xb(3)=w(j,j2,3)
      endif
      if(k2.le.0) then
      xc(1)=x(k,1)
      xc(2)=x(k,2)
      xc(3)=x(k,3)
      else
      xc(1)=w(k,k2,1)
      xc(2)=w(k,k2,2)
      xc(3)=w(k,k2,3)
      endif
      if(l2.le.0) then
      xd(1)=x(l,1)
      xd(2)=x(l,2)
      xd(3)=x(l,3)
      else
      xd(1)=w(l,l2,1)
      xd(2)=w(l,l2,2)
      xd(3)=w(l,l2,3)
      endif
c
c     determine the dihedral angle
c     ----------------------------
c
      xijx=xa(1)-xb(1)
      xkjx=xc(1)-xb(1)
      xklx=xc(1)-xd(1)
      xikx=xijx-xkjx
      xjlx=xklx-xkjx
      xijy=xa(2)-xb(2)
      xkjy=xc(2)-xb(2)
      xkly=xc(2)-xd(2)
      xiky=xijy-xkjy
      xjly=xkly-xkjy
      xijz=xa(3)-xb(3)
      xkjz=xc(3)-xb(3)
      xklz=xc(3)-xd(3)
      xikz=xijz-xkjz
      xjlz=xklz-xkjz
      xmx=xijy*xkjz-xkjy*xijz
      xmy=xijz*xkjx-xkjz*xijx
      xmz=xijx*xkjy-xkjx*xijy
      xnx=xkjy*xklz-xkly*xkjz
      xny=xkjz*xklx-xklz*xkjx
      xnz=xkjx*xkly-xklx*xkjy
      rm2i=1.0d0/(xmx*xmx+xmy*xmy+xmz*xmz)
      rn2i=1.0d0/(xnx*xnx+xny*xny+xnz*xnz)
      rmni=sqrt(rm2i*rn2i)
      cphi=(xmx*xnx+xmy*xny+xmz*xnz)*rmni
      if(cphi.lt.-1.0d0) cphi=-1.0d0
      if(cphi.gt. 1.0d0) cphi= 1.0d0
      phi=acos(cphi)
      s=xkjx*(xmy*xnz-xmz*xny) +xkjy*(xmz*xnx-xmx*xnz)
     + +xkjz*(xmx*xny-xmy*xnx)
      if(s.lt.0.0d0) phi=-phi
c
      argos_diana_torsion=phi
c
      return
      end
