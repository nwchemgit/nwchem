      real*8 function ana_torsion(x,nx,i,j,k,l)
c
c $Id: ana_torsion.F,v 1.1 2000-01-24 23:16:27 d3j191 Exp $
c
      implicit none
c
      integer nx,i,j,k,l
      real*8 x(nx,3)
c
      real*8 xijx,xkjx,xklx,xikx,xjlx
      real*8 xijy,xkjy,xkly,xiky,xjly
      real*8 xijz,xkjz,xklz,xikz,xjlz
      real*8 xmx,xmy,xmz,xnx,xny,xnz,rm2i,rn2i,rmni,cphi,s,phi
c
c     determine the dihedral angle
c     ----------------------------
c
      xijx=x(i,1)-x(j,1)
      xkjx=x(k,1)-x(j,1)
      xklx=x(k,1)-x(l,1)
      xikx=xijx-xkjx
      xjlx=xklx-xkjx
      xijy=x(i,2)-x(j,2)
      xkjy=x(k,2)-x(j,2)
      xkly=x(k,2)-x(l,2)
      xiky=xijy-xkjy
      xjly=xkly-xkjy
      xijz=x(i,3)-x(j,3)
      xkjz=x(k,3)-x(j,3)
      xklz=x(k,3)-x(l,3)
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
      ana_torsion=phi
c
      write(*,'(a)') ' '
      write(*,'(i5,4f12.6)') i,x(i,1),x(i,2),x(i,3)
      write(*,'(i5,4f12.6)') j,x(j,1),x(j,2),x(j,3)
      write(*,'(i5,4f12.6)') k,x(k,1),x(k,2),x(k,3)
      write(*,'(i5,4f12.6)') l,x(l,1),x(l,2),x(l,3),phi
c
      return
      end
