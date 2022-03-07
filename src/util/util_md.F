      subroutine rotate(v,w,g,x,y)
c
c $Id$
c
      implicit none
c
      real*8 v(3),w(3),x(3),y(3),xx(3),t(3,3)
      real*8 small,pi,r,a,b,ca,cb,cg,sa,sb,sg,g
      integer i 
      parameter (small=1.0d-24)
c
c     rotation with angle g around vector from v to w
c     of point x giving result in y
c
      if(abs(w(2)-v(2)).lt.small) then
      if(abs(w(1)-v(1)).lt.small) then
      a=0.0d0
      else
      if(w(1)-v(1).gt.0.0d0) then
      a=2.0d0*datan(1.0d0)
      else
      a=(-2.0d0)*datan(1.0d0)
      endif
      endif
      else
      a=atan(abs(w(1)-v(1))/abs(w(2)-v(2)))
      pi=4.0d0*atan(1.0d0)
      if(w(1)-v(1).gt.0.0d0.and.w(2)-v(2).lt.0.0d0) a=pi-a
      if(w(1)-v(1).lt.0.0d0.and.w(2)-v(2).gt.0.0d0) a=-a
      if(w(1)-v(1).lt.0.0d0.and.w(2)-v(2).lt.0.0d0) a=pi+a
      endif
      r=0.0d0
      do 1 i=1,3
      r=r+(w(i)-v(i))**2
      xx(i)=x(i)-v(i)
    1 continue
      if(r.lt.small) then
      y(1)=x(1)
      y(2)=x(2)
      y(3)=x(3)
      return
      endif
      b=acos((w(3)-v(3))/sqrt(r))
      sa=sin(a)
      ca=cos(a)
      sb=sin(b)
      cb=cos(b)
      sg=sin(g)
      cg=cos(g)
      t(1,1)=ca*ca*cg-sa*ca*cb*sg+sa*ca*cb*sg
     +       +sa*sa*cb*cb*cg+sa*sa*sb*sb
      t(1,2)=(-sa)*ca*cg-ca*ca*cb*sg-sa*sa*cb*sg
     +       +sa*ca*cb*cb*cg+sa*ca*sb*sb
      t(1,3)=ca*sb*sg-sa*sb*cb*cg+sa*sb*cb
      t(2,1)=(-sa)*ca*cg+sa*sa*cb*sg+ca*ca*cb*sg
     +       +sa*ca*cb*cb*cg+sa*ca*sb*sb
      t(2,2)=sa*sa*cg+sa*ca*cb*sg-sa*ca*cb*sg
     +       +ca*ca*cb*cb*cg+ca*ca*sb*sb
      t(2,3)=(-sa)*sb*sg-ca*sb*cb*cg+ca*sb*cb
      t(3,1)=(-ca)*sb*sg-sa*sb*cb*cg+sa*sb*cb
      t(3,2)=sa*sb*sg-ca*sb*cb*cg+ca*sb*cb
      t(3,3)=sb*sb*cg+cb*cb
      do 2 i=1,3
      y(i)=xx(1)*t(i,1)+xx(2)*t(i,2)+xx(3)*t(i,3)+v(i)
    2 continue
      return
      end
      subroutine super(x,nx,mx,y,w,ws,sdev,ny,my,mod,rms0,rms1,
     + xw,mw,nw,ma,na,lpara)
c
c     superimpose x(1:n,1:3) onto y(1:n,1:3)
c
      implicit none
c
#include "msgids.fh"
#include "global.fh"
c
      real*8 zero
      parameter(zero=0.0d0)
c
      integer nx,mx,ny,my,mw,nw,ma,na
      real*8 x(mx,3),y(my,3),w(my),ws(my),sdev(my),xw(mw,ma,3),rms0,rms1
      logical mod,lpara
c
      integer i,j,k,l,nr
      real*8 u(3,3),c(4,4),q(4),b(4),v(4,4),wnorm,wnorms
      real*8 xt(3),yt(3),sd
c
c      write(*,'(a,6f12.6)') 'Super in ',(x(1,i),i=1,3),(x(nx,i),i=1,3)
      wnorm=zero
      wnorms=zero
      do 1 j=1,3
      xt(j)=zero
      yt(j)=zero
    1 continue
      do 2 k=1,nx
      do 3 j=1,3
      xt(j)=xt(j)+w(k)*x(k,j)
    3 continue
    2 continue
      if(lpara) call ga_dgop(mag_d08,xt,3,'+')
      do 4 k=1,ny
      wnorm=wnorm+w(k)
      wnorms=wnorms+ws(k)
      do 5 j=1,3
      yt(j)=yt(j)+w(k)*y(k,j)
    5 continue
    4 continue
c
      rms0=zero
      do 6 j=1,3
      xt(j)=xt(j)/wnorm
      yt(j)=yt(j)/wnorm
      do 7 k=1,nx
      rms0=rms0+ws(k)*(x(k,j)-xt(j)-y(k,j)+yt(j))**2
    7 continue
    6 continue
      if(lpara) call ga_dgop(mag_d09,rms0,1,'+')
      rms0=sqrt(rms0/wnorms)
c
      do 8 i=1,3
      do 9 j=1,3
      u(i,j)=zero
      do 10 k=1,nx
      u(i,j)=u(i,j)+w(k)*(x(k,i)-xt(i))*(y(k,j)-yt(j))
   10 continue
    9 continue
    8 continue
      if(lpara) call ga_dgop(mag_d10,u,9,'+')
c
      c(1,1)=u(1,1)+u(2,2)+u(3,3)
      c(1,2)=u(3,2)-u(2,3)
      c(1,3)=u(1,3)-u(3,1)
      c(1,4)=u(2,1)-u(1,2)
      c(2,2)=u(1,1)-u(2,2)-u(3,3)
      c(2,3)=u(1,2)+u(2,1)
      c(2,4)=u(3,1)+u(1,3)
      c(3,3)=u(2,2)-u(3,3)-u(1,1)
      c(3,4)=u(2,3)+u(3,2)
      c(4,4)=u(3,3)-u(1,1)-u(2,2)
c
      do 11 j=1,3
      do 12 i=j+1,4
      c(i,j)=c(j,i)
   12 continue
   11 continue
c
      call md_jacobi(c,4,4,b,v,nr)
c
      do 13 i=1,4
      q(i)=v(i,4)
   13 continue
c
      u(1,1)=q(1)*q(1)+q(2)*q(2)-q(3)*q(3)-q(4)*q(4)
      u(1,2)=2.0d0*(q(3)*q(2)+q(1)*q(4))
      u(1,3)=2.0d0*(q(4)*q(2)-q(1)*q(3))
      u(2,1)=2.0d0*(q(2)*q(3)-q(1)*q(4))
      u(2,2)=q(1)*q(1)-q(2)*q(2)+q(3)*q(3)-q(4)*q(4)
      u(2,3)=2.0d0*(q(4)*q(3)+q(1)*q(2))
      u(3,1)=2.0d0*(q(2)*q(4)+q(1)*q(3))
      u(3,2)=2.0d0*(q(3)*q(4)-q(1)*q(2))
      u(3,3)=q(1)*q(1)-q(2)*q(2)-q(3)*q(3)+q(4)*q(4)
c
      rms1=zero
      do 14 k=1,nx
      do 15 i=1,3
      q(i)=0.0d0
      do 16 j=1,3
      q(i)=q(i)+u(i,j)*(x(k,j)-xt(j))
   16 continue
   15 continue
      if(mod) then
      do 17 i=1,3
      x(k,i)=q(i)+yt(i)
   17 continue
      endif
      sd=0.0d0
      do 18 j=1,3
      sd=sd+(q(j)-y(k,j)+yt(j))**2
   18 continue
      sdev(k)=sdev(k)+sd
      rms1=rms1+ws(k)*sd
   14 continue
      if(lpara) call ga_dgop(mag_d11,rms1,1,'+')
      rms1=sqrt(rms1/wnorms)
c
      if(mod.and.nw.gt.0) then
      do 19 l=1,nw
      do 20 k=1,na
      do 21 i=1,3
      q(i)=0.0d0
      do 22 j=1,3
      q(i)=q(i)+u(i,j)*(xw(l,k,j)-xt(j))
   22 continue
   21 continue
      do 23 i=1,3
      xw(l,k,i)=q(i)+yt(i)
   23 continue
   20 continue
   19 continue
      endif
c
c      write(*,'(a,6f12.6)') 'Super out',(x(1,i),i=1,3),(x(nx,i),i=1,3)
      return
      end
      subroutine super2(x,ix,nx,mx,y,w,ws,sdev,ny,my,mod,rms0,rms1,
     + xw,mw,nw,ma,na,lpara)
c
c     superimpose x(1:n,1:3) onto y(1:n,1:3)
c
      implicit none
c
#include "msgids.fh"
#include "global.fh"
c
      real*8 zero
      parameter(zero=0.0d0)
c
      integer nx,mx,ny,my,mw,nw,ma,na
      integer ix(mx)
      real*8 x(mx,3),y(my,3),w(my),ws(my),sdev(my),xw(mw,ma,3),rms0,rms1
      logical mod,lpara
c
      integer i,j,k,l,nr
      real*8 u(3,3),c(4,4),q(4),b(4),v(4,4),wnorm,wnorms
      real*8 xt(3),yt(3),sd
c
c      write(*,'(a,6f12.6)') 'Super in ',(x(1,i),i=1,3),(x(nx,i),i=1,3)
      wnorm=zero
      wnorms=zero
      do 1 j=1,3
      xt(j)=zero
      yt(j)=zero
    1 continue
      do 2 k=1,nx
      do 3 j=1,3
      xt(j)=xt(j)+w(ix(k))*x(k,j)
    3 continue
    2 continue
      if(lpara) call ga_dgop(mag_d08,xt,3,'+')
      do 4 k=1,ny
      wnorm=wnorm+w(k)
      wnorms=wnorms+ws(k)
      do 5 j=1,3
      yt(j)=yt(j)+w(k)*y(k,j)
    5 continue
    4 continue
c
      rms0=zero
      do 6 j=1,3
      xt(j)=xt(j)/wnorm
      yt(j)=yt(j)/wnorm
      do 7 k=1,nx
      rms0=rms0+ws(ix(k))*(x(k,j)-xt(j)-y(ix(k),j)+yt(j))**2
    7 continue
    6 continue
      if(lpara) call ga_dgop(mag_d09,rms0,1,'+')
      rms0=sqrt(rms0/wnorms)
c
      do 8 i=1,3
      do 9 j=1,3
      u(i,j)=zero
      do 10 k=1,nx
      u(i,j)=u(i,j)+w(ix(k))*(x(k,i)-xt(i))*(y(ix(k),j)-yt(j))
   10 continue
    9 continue
    8 continue
      if(lpara) call ga_dgop(mag_d10,u,9,'+')
c
      c(1,1)=u(1,1)+u(2,2)+u(3,3)
      c(1,2)=u(3,2)-u(2,3)
      c(1,3)=u(1,3)-u(3,1)
      c(1,4)=u(2,1)-u(1,2)
      c(2,2)=u(1,1)-u(2,2)-u(3,3)
      c(2,3)=u(1,2)+u(2,1)
      c(2,4)=u(3,1)+u(1,3)
      c(3,3)=u(2,2)-u(3,3)-u(1,1)
      c(3,4)=u(2,3)+u(3,2)
      c(4,4)=u(3,3)-u(1,1)-u(2,2)
c
      do 11 j=1,3
      do 12 i=j+1,4
      c(i,j)=c(j,i)
   12 continue
   11 continue
c
      call md_jacobi(c,4,4,b,v,nr)
c
      do 13 i=1,4
      q(i)=v(i,4)
   13 continue
c
      u(1,1)=q(1)*q(1)+q(2)*q(2)-q(3)*q(3)-q(4)*q(4)
      u(1,2)=2.0d0*(q(3)*q(2)+q(1)*q(4))
      u(1,3)=2.0d0*(q(4)*q(2)-q(1)*q(3))
      u(2,1)=2.0d0*(q(2)*q(3)-q(1)*q(4))
      u(2,2)=q(1)*q(1)-q(2)*q(2)+q(3)*q(3)-q(4)*q(4)
      u(2,3)=2.0d0*(q(4)*q(3)+q(1)*q(2))
      u(3,1)=2.0d0*(q(2)*q(4)+q(1)*q(3))
      u(3,2)=2.0d0*(q(3)*q(4)-q(1)*q(2))
      u(3,3)=q(1)*q(1)-q(2)*q(2)-q(3)*q(3)+q(4)*q(4)
c
      rms1=zero
      do 14 k=1,nx
      do 15 i=1,3
      q(i)=0.0d0
      do 16 j=1,3
      q(i)=q(i)+u(i,j)*(x(k,j)-xt(j))
   16 continue
   15 continue
      if(mod) then
      do 17 i=1,3
      x(k,i)=q(i)+yt(i)
   17 continue
      endif
      sd=0.0d0
      do 18 j=1,3
      sd=sd+(q(j)-y(ix(k),j)+yt(j))**2
   18 continue
      sdev(ix(k))=sdev(ix(k))+sd
      rms1=rms1+ws(ix(k))*sd
   14 continue
      if(lpara) call ga_dgop(mag_d11,rms1,1,'+')
      rms1=sqrt(rms1/wnorms)
c
      if(mod.and.nw.gt.0) then
      do 19 l=1,nw
      do 20 k=1,na
      do 21 i=1,3
      q(i)=0.0d0
      do 22 j=1,3
      q(i)=q(i)+u(i,j)*(xw(l,k,j)-xt(j))
   22 continue
   21 continue
      do 23 i=1,3
      xw(l,k,i)=q(i)+yt(i)
   23 continue
   20 continue
   19 continue
      endif
c
c      write(*,'(a,6f12.6)') 'Super out',(x(1,i),i=1,3),(x(nx,i),i=1,3)
      return
      end
      real*8 function angl(x,y,z)
c
      implicit none
c
      real*8 x(3),y(3),z(3)
c
      real*8 xy(3),zy(3),rxy,rzy,phi
      integer i
c
      rxy=0.0d0
      rzy=0.0d0
      do 1 i=1,3
      xy(i)=x(i)-y(i)
      zy(i)=z(i)-y(i)
    1 continue
      rxy=xy(1)*xy(1)+xy(2)*xy(2)+xy(3)*xy(3)
      rzy=zy(1)*zy(1)+zy(2)*zy(2)+zy(3)*zy(3)
      phi=(xy(1)*zy(1)+xy(2)*zy(2)+xy(3)*zy(3))/sqrt(rxy*rzy)
      if(phi.lt.-1.0d0) phi=-1.0d0
      if(phi.gt.1.0d0) phi=1.0d0
      angl=acos(phi)
c
      return
      end
      real*8 function atom_radius(number)
c
      implicit none
      integer number
c
      real*8 radius(0:105)
c
      data radius / 99999.99,
     + 0.35, 1.22, 1.23, 0.89, 0.88, 0.77, 0.70, 0.66, 0.58, 1.60,
     + 1.40, 1.36, 1.25, 1.17, 1.10, 1.04, 0.99, 1.91, 2.03, 1.74,
     + 1.44, 1.32, 1.22, 1.19, 1.17,1.165, 1.16, 1.15, 1.17, 1.25,
     + 1.25, 1.22, 1.21, 1.17, 1.14, 1.98, 2.22, 1.92, 1.62, 1.45,
     + 1.34, 1.29, 1.27, 1.24, 1.25, 1.28, 1.34, 1.41, 1.50, 1.40,
     + 1.41, 1.37, 1.33, 2.09, 2.35, 1.98, 1.69, 1.65, 1.65, 1.64,
     + 1.65, 1.66, 1.65, 1.61, 1.59, 1.59, 1.58, 1.57, 1.56, 1.56,
     + 1.56, 1.44, 1.34, 1.30, 1.28, 1.26, 1.26, 1.29, 1.34, 1.44,
     + 1.55, 1.54, 1.52, 1.53, 1.50, 2.20, 3.24, 2.68, 2.25, 2.16,
     + 1.93, 1.66, 1.57, 1.81, 2.21, 1.43, 1.42, 1.40, 1.39, 1.38,
     + 1.37, 1.36, 1.34, 1.30, 1.30 /
c
      atom_radius=1.1d-01*radius(number)
c
      return
      end
      subroutine povinc(lfn,xmin,xmax,ymin,ymax,zmin,zmax)
c
      implicit none
c
      integer lfn
      real*8 xmin,xmax,ymin,ymax,zmin,zmax
c
      open(unit=lfn,file='camera.inc',form='formatted',status='new',
     + err=9)
      write(lfn,1000)
     + 0.0d0,0.0d0,-1.0d1*max(abs(xmax),abs(ymax),abs(zmax),
     + abs(xmin),abs(ymin),abs(zmin)),
     + 5.0d-1*(xmax+xmin),5.0d-1*(ymax+ymin),5.0d-1*(zmax+zmin)
c     + 0.0d0, 0.0d0, 0.0d0,
c     + 1.3d0, 0.0d0, 0.0d0,
c     + 0.0d0, 0.0d0, 1.0d0,
c     + 0.0d0, 1.0d0, 0.0d0
 1000 format('camera {',/,
     + ' location <',f12.6,',',f12.6,',',f12.6,'>',/,
     + ' look_at <',f12.6,',',f12.6,',',f12.6,'>',/,
     + ' angle 20.0',/,'}')
      write(lfn,1001)  0.0d1, 0.0d1,-5.0d1, 1.0d0,1.0d0,1.0d0
      write(lfn,1001) -1.0d1, 2.0d1,-1.0d1, 1.0d0,1.0d0,1.0d0
      write(lfn,1001)  1.0d1, 2.0d1,-1.0d1, 1.0d0,1.0d0,1.0d0
      write(lfn,1001)  0.0d1, 1.0d1,-2.0d1, 1.0d0,1.0d0,1.0d0
 1001 format('light_source { <',f12.6,',',f12.6,',',f12.6,
     + '> color rgb <',f4.2,',',f4.2,',',f4.2,'> }')
      write(lfn,1002) 0.0d0,0.0d0,0.0d0
 1002 format('background { color rgb <',f4.2,',',f4.2,',',f4.2,'> }')
      close(unit=lfn)
c
    9 continue
      open(unit=lfn,file='colors.inc',form='formatted',status='new',
     + err=99)
      write(lfn,2000)
 2000 format('#declare Colors_Inc_Temp = version ;',/,'#version 2.0 ;')
      write(lfn,2001) 'Gray05',0.05,0.05,0.05
      write(lfn,2001) 'Gray05',0.05,0.05,0.05
      write(lfn,2001) 'Gray10',0.10,0.10,0.10
      write(lfn,2001) 'Gray15',0.15,0.15,0.15
      write(lfn,2001) 'Gray20',0.20,0.20,0.20
      write(lfn,2001) 'Gray25',0.25,0.25,0.25
      write(lfn,2001) 'Gray30',0.30,0.30,0.30
      write(lfn,2001) 'Gray35',0.35,0.35,0.35
      write(lfn,2001) 'Gray40',0.40,0.40,0.40
      write(lfn,2001) 'Gray45',0.45,0.45,0.45
      write(lfn,2001) 'Gray50',0.50,0.50,0.50
      write(lfn,2001) 'Gray55',0.55,0.55,0.55
      write(lfn,2001) 'Gray60',0.60,0.60,0.60
      write(lfn,2001) 'Gray65',0.65,0.65,0.65
      write(lfn,2001) 'Gray70',0.70,0.70,0.70
      write(lfn,2001) 'Gray75',0.75,0.75,0.75
      write(lfn,2001) 'Gray80',0.80,0.80,0.80
      write(lfn,2001) 'Gray85',0.85,0.85,0.85
      write(lfn,2001) 'Gray90',0.90,0.90,0.90
      write(lfn,2001) 'Gray95',0.95,0.95,0.95
      write(lfn,2001) 'DimGray',0.329412,0.329412,0.329412
      write(lfn,2001) 'DimGrey',0.329412,0.329412,0.329412
      write(lfn,2001) 'Gray',0.752941,0.752941,0.752941
      write(lfn,2001) 'Grey',0.752941,0.752941,0.752941
      write(lfn,2001) 'LightGray',0.658824,0.658824,0.658824
      write(lfn,2001) 'LightGrey',0.658824,0.658824,0.658824
      write(lfn,2001) 'VLightGrey',0.80,0.80,0.80
      write(lfn,2001) 'White',1.0,1.0,1.0
      write(lfn,2001) 'Red',1.0,0.0,0.0
      write(lfn,2001) 'Green',0.0,1.0,0.0
      write(lfn,2001) 'Blue',0.0,0.0,1.0
      write(lfn,2001) 'Yellow',1.0,1.0,0.0
      write(lfn,2001) 'Cyan',0.0,1.0,1.0
      write(lfn,2001) 'Magenta',1.0,0.0,1.0
      write(lfn,2001) 'Black',0.0,0.0,0.0
      write(lfn,2001) 'Aquamarine',0.439216,0.858824,0.576471
      write(lfn,2001) 'BlueViolet',0.62352,0.372549,0.623529
      write(lfn,2001) 'Brown',0.647059,0.164706,0.164706
      write(lfn,2001) 'CadetBlue',0.372549,0.623529,0.623529
      write(lfn,2001) 'Coral',1.0,0.498039,0.0
      write(lfn,2001) 'CornflowerBlue',0.258824,0.258824,0.435294
      write(lfn,2001) 'DarkGreen',0.184314,0.309804,0.184314
      write(lfn,2001) 'DarkOliveGreen',0.309804,0.309804,0.184314
      write(lfn,2001) 'DarkOrchid',0.6,0.196078,0.8
      write(lfn,2001) 'DarkSlateBlue',0.419608,0.137255,0.556863
      write(lfn,2001) 'DarkSlateGray',0.184314,0.309804,0.309804
      write(lfn,2001) 'DarkSlateGrey',0.184314,0.309804,0.309804
      write(lfn,2001) 'DarkTurquoise',0.439216,0.576471,0.858824
      write(lfn,2001) 'Firebrick',0.556863,0.137255,0.137255
      write(lfn,2001) 'ForestGreen',0.137255,0.556863,0.137255
      write(lfn,2001) 'Gold',0.8,0.498039,0.196078
      write(lfn,2001) 'Goldenrod',0.858824,0.858824,0.439216
      write(lfn,2001) 'GreenYellow',0.576471,0.858824,0.439216
      write(lfn,2001) 'IndianRed',0.309804,0.184314,0.184314
      write(lfn,2001) 'Khaki',0.623529,0.623529,0.372549
      write(lfn,2001) 'LightBlue',0.74902,0.847059,0.847059
      write(lfn,2001) 'LightSteelBlue',0.560784,0.560784,0.737255
      write(lfn,2001) 'LimeGreen',0.196078,0.8,0.196078
      write(lfn,2001) 'Maroon',0.556863,0.137255,0.419608
      write(lfn,2001) 'MediumAquamarine',0.196078,0.8,0.6
      write(lfn,2001) 'MediumBlue',0.196078,0.196078,0.8
      write(lfn,2001) 'MediumForestGreen',0.419608,0.556863,0.137255
      write(lfn,2001) 'MediumGoldenrod',0.917647,0.917647,0.678431
      write(lfn,2001) 'MediumOrchid',0.576471,0.439216,0.858824
      write(lfn,2001) 'MediumSeaGreen',0.258824,0.435294,0.258824
      write(lfn,2001) 'MediumSlateBlue',0.498039,1.0,0.0
      write(lfn,2001) 'MediumSpringGreen',0.498039,1.0,0.0
      write(lfn,2001) 'MediumTurquoise',0.439216,0.858824,0.858824
      write(lfn,2001) 'MediumVioletRed',0.858824,0.439216,0.576471
      write(lfn,2001) 'MidnightBlue',0.184314,0.184314,0.309804
      write(lfn,2001) 'Navy',0.137255,0.137255,0.556863
      write(lfn,2001) 'NavyBlue',0.137255,0.137255,0.556863
      write(lfn,2001) 'Orange',1,0.5,0.0
      write(lfn,2001) 'OrangeRed',1.0,0.498039,0.0
      write(lfn,2001) 'Orchid',0.858824,0.439216,0.858824
      write(lfn,2001) 'PaleGreen',0.560784,0.737255,0.560784
      write(lfn,2001) 'Pink',0.737255,0.560784,0.560784
      write(lfn,2001) 'Plum',0.917647,0.678431,0.917647
      write(lfn,2001) 'Salmon',0.435294,0.258824,0.258824
      write(lfn,2001) 'SeaGreen',0.137255,0.556863,0.419608
      write(lfn,2001) 'Sienna',0.556863,0.419608,0.137255
      write(lfn,2001) 'SkyBlue',0.196078,0.6,0.8
      write(lfn,2001) 'SlateBlue',0.0,0.498039,1.0
      write(lfn,2001) 'SpringGreen',0.0,1.0,0.498039
      write(lfn,2001) 'SteelBlue',0.137255,0.419608,0.556863
      write(lfn,2001) 'Tan',0.858824,0.576471,0.439216
      write(lfn,2001) 'Thistle',0.847059,0.74902,0.847059
      write(lfn,2001) 'Turquoise',0.678431,0.917647,0.917647
      write(lfn,2001) 'Violet',0.309804,0.184314,0.309804
      write(lfn,2001) 'VioletRed',0.8,0.196078,0.6
      write(lfn,2001) 'Wheat',0.847059,0.847059,0.74902
      write(lfn,2001) 'YellowGreen',0.6,0.8,0.196078
      write(lfn,2001) 'SummerSky',0.22,0.69,0.87
      write(lfn,2001) 'RichBlue',0.35,0.35,0.67
      write(lfn,2001) 'Brass',0.71,0.65,0.26
      write(lfn,2001) 'Copper',0.72,0.45,0.20
      write(lfn,2001) 'Bronze',0.55,0.47,0.14
      write(lfn,2001) 'Bronze2',0.65,0.49,0.24
      write(lfn,2001) 'Silver',0.90,0.91,0.98
      write(lfn,2001) 'BrightGold',0.85,0.85,0.10
      write(lfn,2001) 'OldGold',0.81,0.71,0.23
      write(lfn,2001) 'Feldspar',0.82,0.57,0.46
      write(lfn,2001) 'Quartz',0.85,0.85,0.95
      write(lfn,2001) 'Mica',0.0,0.0,0.0
      write(lfn,2001) 'NeonPink',1.00,0.43,0.78
      write(lfn,2001) 'DarkPurple',0.53,0.12,0.47
      write(lfn,2001) 'NeonBlue',0.30,0.30,1.00
      write(lfn,2001) 'CoolCopper',0.85,0.53,0.10
      write(lfn,2001) 'MandarinOrange',0.89,0.47,0.20
      write(lfn,2001) 'LightWood',0.91,0.76,0.65
      write(lfn,2001) 'MediumWood',0.65,0.50,0.39
      write(lfn,2001) 'DarkWood',0.52,0.37,0.26
      write(lfn,2001) 'SpicyPink',1.00,0.11,0.68
      write(lfn,2001) 'SemiSweetChoc',0.42,0.26,0.15
      write(lfn,2001) 'BakersChoc',0.36,0.20,0.09
      write(lfn,2001) 'Flesh',0.96,0.80,0.69
      write(lfn,2001) 'NewTan',0.92,0.78,0.62
      write(lfn,2001) 'NewMidnightBlue',0.00,0.00,0.61
      write(lfn,2001) 'VeryDarkBrown',0.35,0.16,0.14
      write(lfn,2001) 'DarkBrown',0.36,0.25,0.20
      write(lfn,2001) 'DarkTan',0.59,0.41,0.31
      write(lfn,2001) 'GreenCopper',0.32,0.49,0.46
      write(lfn,2001) 'DkGreenCopper',0.29,0.46,0.43
      write(lfn,2001) 'DustyRose',0.52,0.39,0.39
      write(lfn,2001) 'HuntersGreen',0.13,0.37,0.31
      write(lfn,2001) 'Scarlet',0.55,0.09,0.09
      write(lfn,2002) 'Clear',1.0,1.0,1.0,1.0
 2001 format('#declare ',a,' = color red ',f8.6,' green ',f8.6,
     + ' blue ',f8.6,' ;')
 2002 format('#declare ',a,' = color red ',f8.6,' green ',f8.6,
     + ' blue ',f8.6,' filter ',f8.6,' ;')
      write(lfn,2003)
 2003 format('#declare Plane_Map = 0 ;',/,'#declare Sphere_Map = 1 ;',/,
     + '#declare Cylinder_Map = 2 ;',/,'#declare Torus_Map = 5 ;',/,
     + '#declare Bi   = 2 ;',/,'#declare Norm = 4 ;',/,
     + '#version Colors_Inc_Temp ;')
      close(unit=lfn)
c
   99 continue
      open(unit=lfn,file='plane.inc',form='formatted',status='new',
     + err=999)
      write(lfn,3001)
 3001 format('plane{ <0,1,0>,-1 pigment { White } }') 
      close(unit=lfn)
c
  999 continue
c
      return
      end
      real*8 function torsion(a,b,c,d)
c
      implicit none
c
      real*8 a(3),b(3),c(3),d(3)
c
      torsion=0.0d0
c
      return
      end
      character*255 function atom_color(number)
      implicit none
      integer number
c
      integer indexc(0:105)
c
      data indexc / 12,
     +  7,  5, 14, 12, 13,  0,  1,  2,  6, 12,
     + 13, 15,  9,  6,  8,  3, 13, 12, 16, 16,
     + 12,  9, 12,  9,  9,  8, 12, 10, 10, 10,
     + 12, 12, 12, 12, 10, 12, 12, 12, 12, 12,
     + 12, 12, 12, 12, 12, 12, 16, 12, 12, 12,
     + 12, 12, 11, 12, 12,  8, 12, 12, 12, 12,
     + 12, 12, 12, 12, 12, 12, 12, 12, 12, 12,
     + 12, 12, 12, 12, 12, 12, 12, 16, 17, 16,
     + 12, 12, 12, 12, 12, 12, 12, 12, 12, 12,
     + 12, 12, 12, 12, 12, 12, 12, 12, 12, 12,
     + 12, 12, 12, 12, 12 /
c
      atom_color='white '
      if(indexc(number).eq.0) atom_color='LightGrey '
      if(indexc(number).eq.1) atom_color='SkyBlue '
      if(indexc(number).eq.2) atom_color='Red '
      if(indexc(number).eq.3) atom_color='Yellow '
      if(indexc(number).eq.4) atom_color='White '
      if(indexc(number).eq.5) atom_color='Pink '
      if(indexc(number).eq.6) atom_color='Goldenrod '
      if(indexc(number).eq.7) atom_color='Blue '
      if(indexc(number).eq.8) atom_color='Orange '
      if(indexc(number).eq.9) atom_color='Gray85 '
      if(indexc(number).eq.10) atom_color='Brown '
      if(indexc(number).eq.11) atom_color='Purple '
      if(indexc(number).eq.12) atom_color='SpicyPink '
      if(indexc(number).eq.13) atom_color='Green '
      if(indexc(number).eq.14) atom_color='Firebrick '
      if(indexc(number).eq.15) atom_color='DarkGreen '
      if(indexc(number).eq.16) atom_color='Silver '
      if(indexc(number).eq.17) atom_color='Gold '
c
      return
      end
      subroutine md_jacobi(a,n,na,d,v,nrot)
c
c     compute eigenvectors and eigenvalues for real symmetric 
c     matrix using the Jacobi diagonalization
c
      implicit none
c
      integer nrmax
      parameter(nrmax=100)
c
      real*8 zero,half,one,two
      parameter(zero=0.0d0)
      parameter(half=0.5d0)
      parameter(one=1.0d0)
      parameter(two=2.0d0)
c
      integer n,na,nrot
      real*8 a(na,na),d(na),v(na,na)
      real*8 at,b,dma,q
c
      integer i,j,k,l
      real*8 c,s,t,sum,temp
c
      do 1 i=1,n
      do 2 j=1,n
      v(i,j)=zero
    2 continue
      v(i,i)=one
      d(i)=a(i,i)
    1 continue
c
      nrot=0
      do 3 l=1,nrmax
      nrot=nrot+1
      sum=zero
      do 4 i=1,n-1
      do 5 j=i+1,n
      sum=sum+abs(a(i,j))
    5 continue
    4 continue
      if(sum.eq.zero) then
      do 6 i=1,n-1
      do 7 j=i+1,n
      if(d(i).gt.d(j)) then
      temp=d(i)
      d(i)=d(j)
      d(j)=temp
      do 8 k=1,n
      temp=v(k,i)
      v(k,i)=v(k,j)
      v(k,j)=temp
    8 continue
      endif
    7 continue
    6 continue
      return
      endif
      do 9 j=2,n
      do 10 i=1,j-1
      b=a(i,j)
      if(abs(b).gt.zero) then
      dma=d(j)-d(i)
      if(abs(dma)+abs(b).le.abs(dma)) then
      t=b/dma
      else
      q=half*dma/b
      t=sign(one/(abs(q)+sqrt(one+q*q)),q)
      endif
      c=one/sqrt(t*t+one)
      s=t*c
      a(i,j)=zero
      do 11 k=1,i-1
      at=c*a(k,i)-s*a(k,j)
      a(k,j)=s*a(k,i)+c*a(k,j)
      a(k,i)=at
   11 continue
      do 12 k=i+1,j-1
      at=c*a(i,k)-s*a(k,j)
      a(k,j)=s*a(i,k)+c*a(k,j)
      a(i,k)=at
   12 continue
      do 13 k=j+1,n
      at=c*a(i,k)-s*a(j,k)
      a(j,k)=s*a(i,k)+c*a(j,k)
      a(i,k)=at
   13 continue
      do 14 k=1,n
      at=c*v(k,i)-s*v(k,j)
      v(k,j)=s*v(k,i)+c*v(k,j)
      v(k,i)=at
   14 continue
      at=c*c*d(i)+s*s*d(j)-two*c*s*b
      d(j)=s*s*d(i)+c*c*d(j)+two*c*s*b
      d(i)=at
      endif
   10 continue
    9 continue
    3 continue
c
      call md_abort('md_jacobi: maximum iterations reached',0)
c
      return
      end
      subroutine swatch(today,now)
c
      implicit none
c
      character*10 today,now
c
#if defined(LINUX)
      character*26 string
#endif
#if defined(IBM)
      character*26 string
#endif
#if defined(SP1) || defined(SOLARIS)
      character*26 string
#endif
c
      today='00/00/00  '
      now='00:00:00'
c
#if defined(IBM)
      call fdate(string)
      if(string(4:6).eq.'Jan') today(1:2)='01'
      if(string(4:6).eq.'Feb') today(1:2)='02'
      if(string(4:6).eq.'Mar') today(1:2)='03'
      if(string(4:6).eq.'Apr') today(1:2)='04'
      if(string(4:6).eq.'May') today(1:2)='05'
      if(string(4:6).eq.'Jun') today(1:2)='06'
      if(string(4:6).eq.'Jul') today(1:2)='07'
      if(string(4:6).eq.'Aug') today(1:2)='08'
      if(string(4:6).eq.'Sep') today(1:2)='09'
      if(string(4:6).eq.'Oct') today(1:2)='10'
      if(string(4:6).eq.'Nov') today(1:2)='11'
      if(string(4:6).eq.'Dec') today(1:2)='12'
      today(7:8)=string(8:9)
      today(4:5)=string(1:2)
      now=string(11:20)
#endif
#if defined(SP1) || defined(SOLARIS)
      call util_date(string)
      if(string(5:7).eq.'Jan') today(1:2)='01'
      if(string(5:7).eq.'Feb') today(1:2)='02'
      if(string(5:7).eq.'Mar') today(1:2)='03'
      if(string(5:7).eq.'Apr') today(1:2)='04'
      if(string(5:7).eq.'May') today(1:2)='05'
      if(string(5:7).eq.'Jun') today(1:2)='06'
      if(string(5:7).eq.'Jul') today(1:2)='07'
      if(string(5:7).eq.'Aug') today(1:2)='08'
      if(string(5:7).eq.'Sep') today(1:2)='09'
      if(string(5:7).eq.'Oct') today(1:2)='10'
      if(string(5:7).eq.'Nov') today(1:2)='11'
      if(string(5:7).eq.'Dec') today(1:2)='12'
      today(7:8)=string(23:24)
      today(4:5)=string(9:10)
      now=string(11:20)
#endif
#if defined(LINUX)
      call util_date(string)
      if(string(5:7).eq.'Jan') today(1:2)='01'
      if(string(5:7).eq.'Feb') today(1:2)='02'
      if(string(5:7).eq.'Mar') today(1:2)='03'
      if(string(5:7).eq.'Apr') today(1:2)='04'
      if(string(5:7).eq.'May') today(1:2)='05'
      if(string(5:7).eq.'Jun') today(1:2)='06'
      if(string(5:7).eq.'Jul') today(1:2)='07'
      if(string(5:7).eq.'Aug') today(1:2)='08'
      if(string(5:7).eq.'Sep') today(1:2)='09'
      if(string(5:7).eq.'Oct') today(1:2)='10'
      if(string(5:7).eq.'Nov') today(1:2)='11'
      if(string(5:7).eq.'Dec') today(1:2)='12'
      today(7:8)=string(23:24)
      today(4:5)=string(9:10)
      now=string(11:20)
#endif
      if(today(4:4).eq.' ') today(4:4)='0'
      return
      end
      subroutine matinv(a,n,ndim)
c
      implicit none

      integer maxdim
      real*8 zero,small,one
      parameter(maxdim=3)
      parameter(zero=0.0d0)
      parameter(small=1.0d-6)
      parameter(one=1.0d0)
c
      integer n,ndim
      real*8 a(ndim,ndim)
      integer ia(2,maxdim),ib(maxdim),ic(maxdim)
      real*8 d(maxdim)
c
      integer idim,i,j,k,l,m
      real*8 b,e
c
      if(ndim.gt.maxdim) call md_abort('matinv dimension error',0)
c
      do 1 idim=1,n
      ia(1,idim)=0
      ia(2,idim)=0
    1 continue
c
      do 9 idim=1,n
      b=zero
      do 3 l=1,n
      do 4 m=1,n
      if(ia(1,l).ne.1.and.ia(2,m).ne.1) then
      e=dabs(a(l,m))
      if(e.ge.b) then
      i=l
      k=m
      endif
    8 b=dmax1(b,e)
      endif
    4 continue
    3 continue
      ia(1,i)=1
      ia(2,k)=1
      ib(k)=i
      ic(i)=k
      b=a(i,k)
c
      if(dabs(b).lt.small) call md_abort('arg_matinv singular matrix',0)
      a(i,k)=one/b
      do 6 l=1,n
      if(l.ne.k) a(i,l)=-a(i,l)/b
  6   continue
      do 5 l=1,n
      do 7 m=1,n
      if(l.ne.i.and.m.ne.k) a(l,m)=a(l,m)+a(l,k)*a(i,m)
    7 continue
  5   continue
      do 11 l=1,n
      if(l.ne.i) a(l,k)=a(l,k)/b
  11  continue
  9   continue
c
      do 15 l=1,n
      do 13 j=1,n
      k=ib(j)
      d(j)=a(k,l)
   13 continue
      do 14 j=1,n
      a(j,l)=d(j)
   14 continue
  15  continue
c
      do 16 l=1,n
      do 17 j=1,n
      k=ic(j)
      d(j)=a(l,k)
   17 continue
      do 18 j=1,n
      a(l,j)=d(j)
   18 continue
  16  continue
c
      return
      end
      logical function frequency(istep,nstep)
c
      implicit none
c
      integer istep,nstep
c
      if(nstep.le.0) then
      frequency=.false.
      else
      frequency=mod(istep,nstep).eq.0
      endif
c
      return
      end
      subroutine rolex(elaps,cputim)
c
      implicit none
c
      real*8 elaps,cputim
c
      real*8 util_wallsec,util_cpusec
      external util_wallsec,util_cpusec
c
      elaps=util_wallsec()
      cputim=util_cpusec()
c
      return
      end
      subroutine timer_init()
c
      implicit none
c
      call timer(0,0)
c
      return
      end
      subroutine timer_reset(itime)
c
      implicit none
c
      integer itime
c
      call timer(itime,-2)
c
      return
      end
      subroutine timer_start(itime)
c
      implicit none
c
      integer itime
c
      call timer(itime,0)
c
      return
      end
      subroutine timer_stop(itime)
c
      implicit none
c
      integer itime
c
      call timer(itime,1)
c
      return
      end
      subroutine timer(itime,iopt)
c
      implicit none
c
      integer itime,iopt
c
      real*8 elaps,cputim
c
      integer mtime
      parameter(mtime=250)
      integer ncall(250)
      real*8 ttime(250,3),ctime(250,3)
      common/tim/ncall,ttime,ctime
c
      integer i
c
      call rolex(elaps,cputim)
c
      if(itime.eq.0) then
      do 1 i=1,mtime
      ncall(i)=0
      ctime(i,1)=0.0d0
      ttime(i,1)=0.0d0
      ctime(i,2)=0.0d0
      ttime(i,2)=0.0d0
      ctime(i,3)=1.0d9
      ttime(i,3)=1.0d9
    1 continue
      elseif(itime.le.0.or.itime.gt.mtime) then
      call md_abort('Timer index out of range',0)
      elseif(iopt.eq.-2) then
      ncall(itime)=0
      ctime(itime,1)=0.0d0
      ttime(itime,1)=0.0d0
      ctime(itime,2)=0.0d0
      ttime(itime,2)=0.0d0
      ctime(itime,3)=1.0d9
      ttime(itime,3)=1.0d9
      elseif(iopt.eq.-1) then
      ncall(itime)=0
      ctime(itime,1)=-cputim
      ttime(itime,1)=-elaps
      ctime(itime,2)=-cputim
      ttime(itime,2)=-elaps
      ctime(itime,3)=1.0d9
      ttime(itime,3)=1.0d9
      elseif(iopt.eq.0) then
      ctime(itime,1)=ctime(itime,1)-cputim
      ttime(itime,1)=ttime(itime,1)-elaps
      ctime(itime,2)=-cputim
      ttime(itime,2)=-elaps
      elseif(iopt.eq.1) then
      ncall(itime)=ncall(itime)+1
      ctime(itime,1)=ctime(itime,1)+cputim
      ttime(itime,1)=ttime(itime,1)+elaps
      ctime(itime,2)=ctime(itime,2)+cputim
      ttime(itime,2)=ttime(itime,2)+elaps
      ctime(itime,3)=min(ctime(itime,2),ctime(itime,3))
      ttime(itime,3)=min(ttime(itime,2),ttime(itime,3))
      else
      call md_abort('Unimplemented timer option',0)
      endif
c
      return
      end
      real*8 function timer_cpu(itime)
c
      implicit none
c
      integer itime
c
      integer mtime
      parameter(mtime=250)
      integer ncall(250)
      real*8 ttime(250,3),ctime(250,3)
      common/tim/ncall,ttime,ctime
c
      if(itime.le.0.or.itime.gt.mtime)
     + call md_abort('Illegal timer index',0)
c
      if(ncall(itime).le.0) then
      timer_cpu=0.0d0
      else
      timer_cpu=ctime(itime,2)
      endif
c
      return
      end
      real*8 function timer_wall(itime)
c
      implicit none
c
      integer itime
c
      integer mtime
      parameter(mtime=250)
      integer ncall(250)
      real*8 ttime(250,3),ctime(250,3)
      common/tim/ncall,ttime,ctime
c
      if(itime.le.0.or.itime.gt.mtime)
     + call md_abort('Illegal timer index',0)
c
      if(ncall(itime).le.0) then
      timer_wall=0.0d0
      else
      timer_wall=ttime(itime,2)
      endif
c
      return
      end
      integer function timer_calls(itime)
c
      integer itime
c
      integer mtime
      parameter(mtime=250)
      integer ncall(250)
      real*8 ttime(250,3),ctime(250,3)
      common/tim/ncall,ttime,ctime
c
      if(itime.le.0.or.itime.gt.mtime)
     + call md_abort('Illegal timer index',0)
c
      timer_calls=ncall(itime)
c
      return
      end
      real*8 function timer_cpu_minimum(itime)
c
      implicit none
c
      integer itime
c
      integer mtime
      parameter(mtime=250)
      integer ncall(250)
      real*8 ttime(250,3),ctime(250,3)
      common/tim/ncall,ttime,ctime
c
      if(itime.le.0.or.itime.gt.mtime)
     + call md_abort('Illegal timer index',0)
c
      if(ncall(itime).le.0) then
      timer_cpu_minimum=0.0d0
      else
      timer_cpu_minimum=ctime(itime,3)
      endif
c
      return
      end
      real*8 function timer_wall_minimum(itime)
c
      implicit none
c
      integer itime
c
      integer mtime
      parameter(mtime=250)
      integer ncall(250)
      real*8 ttime(250,3),ctime(250,3)
      common/tim/ncall,ttime,ctime
c
      if(itime.le.0.or.itime.gt.mtime)
     + call md_abort('Illegal timer index',0)
c
      if(ncall(itime).le.0) then
      timer_wall_minimum=0.0d0
      else
      timer_wall_minimum=ttime(itime,3)
      endif
c
      return
      end
      real*8 function timer_cpu_average(itime)
c
      implicit none
c
      integer itime
c
      integer mtime
      parameter(mtime=250)
      integer ncall(250)
      real*8 ttime(250,3),ctime(250,3)
      common/tim/ncall,ttime,ctime
c
      if(itime.le.0.or.itime.gt.mtime)
     + call md_abort('Illegal timer index',0)
c
      if(ncall(itime).le.0) then
      timer_cpu_average=0.0d0
      else
      timer_cpu_average=ctime(itime,1)/dble(ncall(itime))
      endif
c
      return
      end
      real*8 function timer_wall_average(itime)
c
      implicit none
c
      integer itime
c
      integer mtime
      parameter(mtime=250)
      integer ncall(250)
      real*8 ttime(250,3),ctime(250,3)
      common/tim/ncall,ttime,ctime
c
      if(itime.le.0.or.itime.gt.mtime)
     + call md_abort('Illegal timer index',0)
c
      if(ncall(itime).le.0) then
      timer_wall_average=0.0d0
      else
      timer_wall_average=ttime(itime,1)/dble(ncall(itime))
      endif
c
      return
      end
      real*8 function timer_cpu_total(itime)
c
      implicit none
c
      integer itime
c
      integer mtime
      parameter(mtime=250)
      integer ncall(250)
      real*8 ttime(250,3),ctime(250,3)
      common/tim/ncall,ttime,ctime
c
      if(itime.le.0.or.itime.gt.mtime)
     + call md_abort('Illegal timer index',0)
c
      if(ncall(itime).le.0) then
      timer_cpu_total=0.0d0
      else
      timer_cpu_total=ctime(itime,1)
      endif
c
      return
      end
      real*8 function timer_wall_total(itime)
c
      implicit none
c
      integer itime
c
      integer mtime
      parameter(mtime=250)
      integer ncall(250)
      real*8 ttime(250,3),ctime(250,3)
      common/tim/ncall,ttime,ctime
c
      if(itime.le.0.or.itime.gt.mtime)
     + call md_abort('Illegal timer index',0)
c
      if(ncall(itime).le.0) then
      timer_wall_total=0.0d0
      else
      timer_wall_total=ttime(itime,1)
      endif
c
      return
      end
      subroutine error(lauto,lapprox,maxacf,data,ndata,
     + aver,drift,stderr,corerr,ratio)
c
      implicit none
c
      real*8 zero,one
      parameter(zero=0.0d0)
      parameter(one=1.0d0)
c
#include "mafdecls.fh"
c
      logical lauto,lapprox
      integer ndata,lenacf,maxacf
      real*8 data(ndata),aver,drift,stderr,corerr,ratio
c
      integer i,i_acf,l_acf
      real*8 dsum,ddsum,dtsum,tsum,ttsum,dstep
c
c      integer nacf,kapprx(15),iapprx,klarge
c      real*8 warg
c      real*8 data(ndata),acf(ndata),approx(15),cdac(15),weight
c
c      integer i,j,k,l,m,nacfa,nfunc
c      real*8 dsum,ddsum,dtsum,tsum,ttsum,dstep,dfsum,dvar
c      real*8 cdap(15,15),cdaq(15,15),cdad(15),rfact
c      real*8 xappm,xappi,xappj,wterm,wsum1,wsum2,cdawgt
c
c     arg_error iopt   = 0 : average
c                            drift
c                            standard error
c                      = 1 : average
c                            drift
c                            standard error
c                            autocorrelation function
c                            correlation error from actual acf
c                            sampling ration from actual acf
c                      = 2 : average
c                            drift
c                            standard error
c                            autocorrelation function
c                            correlation error from approximated acf
c                            sampling ration from approximated acf
c               data   : 1-dimensional array with data
c               ndata  : number of data
c               nacf   : length of autocorrelation function
c
c               aver   : average
c               stderr : standard error
c               drift  : drift
c               acf    : autocorrelation function
c               corerr : correlated error
c               ratio  : sampling ratio
c
      lenacf=maxacf
c
      dsum=zero
      ddsum=zero
      dtsum=zero
      tsum=zero
      ttsum=zero
      do 1 i=1,ndata
      dstep=dble(i)
      dsum=dsum+data(i)
      ddsum=ddsum+data(i)*data(i)
      dtsum=dtsum+dstep*data(i)
      tsum=tsum+dstep
      ttsum=ttsum+dstep*dstep
    1 continue
c
c     average, drift and standard error
c
      aver=dsum/dble(ndata)
      drift=(dble(ndata)*dtsum-dsum*tsum)/(dble(ndata)*ttsum-tsum*tsum)
      stderr=sqrt(abs((ddsum/dble(ndata)-dsum*dsum/dble(ndata*ndata)))/
     + (dble(ndata-1)))
c
      corerr=stderr
      ratio=one
c
      if(.not.lauto) return
c
      if(.not.ma_push_get(mt_dbl,lenacf,'acf',l_acf,i_acf))
     + call md_abort('Failed to allocate memory for acf',0)
c
      call auto_corr(data,ndata,aver,dbl_mb(i_acf),lenacf,ratio)
c
      if(.not.ma_pop_stack(l_acf))
     + call md_abort('Failed to deallocate memory for acf',0)
c
      corerr=ratio*stderr
c
      return
      end
      subroutine auto_corr(data,ndata,aver,acf,lacf,ratio)
c
      implicit none
c
      real*8 zero,half,one,two
      parameter(zero=0.0d0)
      parameter(half=5.0d-1)
      parameter(one=1.0d0)
      parameter(two=2.0d0)
c
      integer ndata,lacf,nacf
      real*8 data(ndata),acf(lacf),aver,ratio
c
      integer i,j
      real*8 dsum,dvar
c
      nacf=min(ndata,lacf)
c
      dsum=zero
      do 1 i=1,ndata
      dsum=dsum+(data(i)-aver)**2
    1 continue
      dvar=dble(ndata)/dsum
      do 2 i=1,nacf-1
      dsum=zero
      do 3 j=1,ndata-i
      dsum=dsum+(data(j)-aver)*(data(i+j)-aver)
    3 continue
      acf(i)=dvar*(dsum/dble(ndata-i))*(one-dble(i-1)/dble(ndata-2))
    2 continue
c
      ratio=half
      do 4 i=1,nacf-1
      ratio=ratio+(one-(dble(i)/dble(ndata))**2)*abs(acf(i))
    4 continue
      ratio=sqrt(two*abs(ratio))
      if(ratio.lt.one) ratio=one
c
      lacf=nacf
c
      return
      end
      subroutine acf_approx(acf,acfapp,lacf)
c
      implicit none
c
      integer lacf
      real*8 acf(lacf),acfapp(lacf)
c
c      integer kapprox(15)
c
c      rfact=two*sqrt(dble(klarge))/dble(nacfa-1)
c      wsum1=zero
c      wsum2=zero
c      warg=zero
c      do 5 i=1,nacfa-1
c      acf(i)=abs(acf(i))
c      if(acf(i).gt.zero) then
c      wterm=exp(weight*dble(nacfa-i)/dble(nacfa-1))
c      wsum1=wsum1+wterm
c      wsum2=wsum2+wterm*log(acf(i))/(rfact*dble(i))
c      endif
c    5 continue
c      if(abs(wsum1).gt.small) warg=wsum2/wsum1
c      do 55 i=1,nacfa-1
c      acf(i)=acf(i)-exp(warg*dble(i)*rfact)
c   55 continue
c      nfunc=iapprx
c      if(nfunc.le.0) nfunc=15
c      do 6 k=1,nfunc
c      do 7 l=1,nfunc
c      cdap(k,l)=zero
c      do 8 m=1,nacfa-1
c      xappm=dble(m)*rfact
c      cdawgt=exp(weight*dble(nacfa-m)/dble(nacfa-1))
c      cdap(k,l)=cdap(k,l)+cdawgt*exp((-xappm)*xappm)*
c     + approx(k)*approx(l)*(xappm**(kapprx(k)+kapprx(l)))
c    8 continue
c      cdaq(k,l)=cdap(k,l)
c    7 continue
c    6 continue
c      do 9 i=1,nfunc
c      cdad(i)=zero
c      do 10 j=1,nacfa-1
c      xappj=dble(j)*rfact
c      cdawgt=exp(weight*dble(nacfa-j)/dble(nacfa-1))
c      cdad(i)=cdad(i)+cdawgt*exp((-half)*xappj*xappj)*
c     + acf(j)*approx(i)*(xappj**kapprx(i))
c   10 continue
c    9 continue
c      do 11 k=1,nfunc
c      do 12 i=k,nfunc
c      cdaq(i,k)=cdap(i,k)
c      do 13 l=1,k-1
c      cdaq(i,k)=cdaq(i,k)-cdaq(i,l)*cdaq(l,k)
c   13 continue
c   12 continue
c      do 14 i=k+1,nfunc
c      cdaq(k,i)=cdap(k,i)
c      do 15 l=1,k-1
c      cdaq(k,i)=cdaq(k,i)-cdaq(k,l)*cdaq(l,i)
c   15 continue
c      cdaq(k,i)=cdaq(k,i)/cdaq(k,k)
c   14 continue
c   11 continue
c      do 16 j=1,nfunc
c      cdac(j)=cdad(j)
c      do 17 i=1,j-1
c      cdac(j)=cdac(j)-cdaq(j,i)*cdac(i)
c   17 continue
c      cdac(j)=cdac(j)/cdaq(j,j)
c   16 continue
c      do 18 k=1,nfunc
c      j=nfunc+1-k
c      do 19 i=j+1,nfunc
c      cdac(j)=cdac(j)-cdaq(j,i)*cdac(i)
c   19 continue
c   18 continue
c      do 20 i=1,nacfa-1
c      xappi=dble(i)*rfact
c      acf(i)=exp(warg*xappi)
c      do 21 j=1,nfunc
c      acf(i)=acf(i)+exp((-half)*xappi*xappi)*
c     + approx(j)*cdac(j)*(xappi**kapprx(j))
c   21 continue
c   20 continue
c
c     autocorrelation function upto lag min(nacf,ndata)
c
c      if(iopt.gt.0) then
c      warg=zero
c      nacfa=nacf
c      if(nacfa.gt.ndata) nacfa=ndata
c      dfsum=zero
c      do 2 i=1,ndata
c      dfsum=dfsum+(data(i)-aver)**2
c    2 continue
c      dvar=dble(ndata)/dfsum
c      do 3 i=1,nacfa-1
c      dfsum=zero
c      do 4 j=1,ndata-i
c      dfsum=dfsum+(data(j)-aver)*(data(i+j)-aver)
c    4 continue
c      acf(i)=dvar*(dfsum/dble(ndata-i))*(one-dble(i-1)/dble(ndata-2))
c    3 continue
c      endif
cc
cc     approximate acf here
cc
c      if(iopt.gt.1) then
c      rfact=two*sqrt(dble(klarge))/dble(nacfa-1)
c      wsum1=zero
c      wsum2=zero
c      warg=zero
c      do 5 i=1,nacfa-1
c      acf(i)=abs(acf(i))
c      if(acf(i).gt.zero) then
c      wterm=exp(weight*dble(nacfa-i)/dble(nacfa-1))
c      wsum1=wsum1+wterm
c      wsum2=wsum2+wterm*log(acf(i))/(rfact*dble(i))
c      endif
c    5 continue
c      if(abs(wsum1).gt.small) warg=wsum2/wsum1
c      do 55 i=1,nacfa-1
c      acf(i)=acf(i)-exp(warg*dble(i)*rfact)
c   55 continue
c      nfunc=iapprx
c      if(nfunc.le.0) nfunc=15
c      do 6 k=1,nfunc
c      do 7 l=1,nfunc
c      cdap(k,l)=zero
c      do 8 m=1,nacfa-1
c      xappm=dble(m)*rfact
c      cdawgt=exp(weight*dble(nacfa-m)/dble(nacfa-1))
c      cdap(k,l)=cdap(k,l)+cdawgt*exp((-xappm)*xappm)*
c     + approx(k)*approx(l)*(xappm**(kapprx(k)+kapprx(l)))
c    8 continue
c      cdaq(k,l)=cdap(k,l)
c    7 continue
c    6 continue
c      do 9 i=1,nfunc
c      cdad(i)=zero
c      do 10 j=1,nacfa-1
c      xappj=dble(j)*rfact
c      cdawgt=exp(weight*dble(nacfa-j)/dble(nacfa-1))
c      cdad(i)=cdad(i)+cdawgt*exp((-half)*xappj*xappj)*
c     + acf(j)*approx(i)*(xappj**kapprx(i))
c   10 continue
c    9 continue
c      do 11 k=1,nfunc
c      do 12 i=k,nfunc
c      cdaq(i,k)=cdap(i,k)
c      do 13 l=1,k-1
c      cdaq(i,k)=cdaq(i,k)-cdaq(i,l)*cdaq(l,k)
c   13 continue
c   12 continue
c      do 14 i=k+1,nfunc
c      cdaq(k,i)=cdap(k,i)
c      do 15 l=1,k-1
c      cdaq(k,i)=cdaq(k,i)-cdaq(k,l)*cdaq(l,i)
c   15 continue
c      cdaq(k,i)=cdaq(k,i)/cdaq(k,k)
c   14 continue
c   11 continue
c      do 16 j=1,nfunc
c      cdac(j)=cdad(j)
c      do 17 i=1,j-1
c      cdac(j)=cdac(j)-cdaq(j,i)*cdac(i)
c   17 continue
c      cdac(j)=cdac(j)/cdaq(j,j)
c   16 continue
c      do 18 k=1,nfunc
c      j=nfunc+1-k
c      do 19 i=j+1,nfunc
c      cdac(j)=cdac(j)-cdaq(j,i)*cdac(i)
c   19 continue
c   18 continue
c      do 20 i=1,nacfa-1
c      xappi=dble(i)*rfact
c      acf(i)=exp(warg*xappi)
c      do 21 j=1,nfunc
c      acf(i)=acf(i)+exp((-half)*xappi*xappi)*
c     + approx(j)*cdac(j)*(xappi**kapprx(j))
c   21 continue
c   20 continue
c      endif
cc
cc     sampling ratio
cc
c      ratio=one
c      if(iopt.gt.0) then
c      ratio=half
c      do 22 i=1,nacfa-1
c      ratio=ratio+(one-(dble(i)/dble(ndata))**2)*abs(acf(i))
c   22 continue
c      ratio=sqrt(two*abs(ratio))
c      if(ratio.lt.one) ratio=one
c      endif
cc
cc
c      corerr=ratio*stderr
cc
      return
      end
      logical function md_zmat(x1,x2,x3,x4,d,a,t)
c
      implicit none
c
      real*8 x1(3),x2(3),x3(3),x4(3),d,a,t
c
      real*8 p(3),q(3),s(3),v(3),w(3)
      real*8 r
      integer i
c
      do 1 i=1,3
      p(i)=x3(i)-x2(i)
      q(i)=x4(i)-x2(i)
    1 continue
      v(1)=p(2)*q(3)-p(3)*q(2)+x2(1)
      v(2)=p(3)*q(1)-p(1)*q(3)+x2(2)
      v(3)=p(1)*q(2)-p(2)*q(1)+x2(3)
c
      r=sqrt(p(1)*p(1)+p(2)*p(2)+p(3)*p(3))
      do 2 i=1,3
      s(i)=d*p(i)/r+x2(i)
    2 continue
c
      call rotate(x2,v,a,s,w)
      call rotate(x3,x2,t,w,s)
c
      do 3 i=1,3
      x1(i)=s(i)
    3 continue
c
      md_zmat=.true.
      return
      end
      subroutine md_abort(string, icode)
      implicit none
#include "global.fh"
#include "stdio.fh"
      character*(*) string
      character*255 card
      integer icode
      if(ga_nodeid().eq.0) then
      write(luout,1000) 0,string,icode
 1000 format(/,1x,10('*'),/,' * ',i3,': ',a,i5,/,1x,10('*'))
      card=' '
      else
      write(card,1001) ga_nodeid(),string,icode
 1001 format(' * ',i3,': ',a,i5,' * ')
      endif
      call ga_error(card,icode)
      return
      end
