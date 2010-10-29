      real*8 function argos_diana_angle(x,nx,w,nw,na,i,i2,j,j2,k,k2)
c
c $Id$
c
      implicit none
c
      integer nx,nw,na,i,i2,j,j2,k,k2
      real*8 x(nx,3),w(nw,na,3)
c
      real*8 one
      parameter(one=1.0d0)
c
      real*8 xsijx,xsijy,xsijz,xskjx,xskjy,xskjz
      real*8 rsij2,rskj2,rsij2i,rskj2i,rsikji,cphi
c
      real*8 xa(3),xb(3),xc(3)
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
c
c     determine the angle
c     -------------------
c
      xsijx=xa(1)-xb(1)
      xskjx=xc(1)-xb(1)
      xsijy=xa(2)-xb(2)
      xskjy=xc(2)-xb(2)
      xsijz=xa(3)-xb(3)
      xskjz=xc(3)-xb(3)
c
      rsij2=xsijx*xsijx+xsijy*xsijy+xsijz*xsijz
      rskj2=xskjx*xskjx+xskjy*xskjy+xskjz*xskjz
      cphi=xsijx*xskjx+xsijy*xskjy+xsijz*xskjz
      rsij2i=one/rsij2
      rskj2i=one/rskj2
      rsikji=one/sqrt(rsij2*rskj2)
      cphi=cphi*rsikji
      if(cphi.lt.-one) cphi=-one
      if(cphi.gt. one) cphi= one
      argos_diana_angle=acos(cphi)
c
      return
      end
