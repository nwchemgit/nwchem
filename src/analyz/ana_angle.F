      real*8 function ana_angle(x,nx,i,j,k)
c
c $Id: ana_angle.F,v 1.1 2000-03-16 20:27:40 d3j191 Exp $
c
      implicit none
c
      integer nx,i,j,k
      real*8 x(nx,3)
c
      real*8 one
      parameter(one=1.0d0)
c
      real*8 xsijx,xsijy,xsijz,xskjx,xskjy,xskjz
      real*8 rsij2,rskj2,rsij2i,rskj2i,rsikji,cphi
c
c     determine the angle
c     -------------------
c
      xsijx=x(i,1)-x(j,1)
      xskjx=x(k,1)-x(j,1)
      xsijy=x(i,2)-x(j,2)
      xskjy=x(k,2)-x(j,2)
      xsijz=x(i,3)-x(j,3)
      xskjz=x(k,3)-x(j,3)
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
      ana_angle=acos(cphi)
c
      return
      end
