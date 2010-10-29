      real*8 function argos_diana_bond(x,nx,w,nw,na,i,i2,j,j2,box)
c
c $Id$
c
      implicit none
c
      integer nx,nw,na,i,i2,j,j2
      real*8 x(nx,3),w(nw,na,3),box(3)
c
      real*8 xa(3),xb(3),dx(3)
      integer ix
c
c     determine the bond length
c     -------------------------
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
c
      do 1 ix=1,3
      dx(ix)=xa(ix)-xb(ix)
      if(dx(ix).gt.box(ix)) dx(ix)=dx(ix)-box(ix)
      if(dx(ix).lt.-box(ix)) dx(ix)=dx(ix)+box(ix)
    1 continue
c      argos_diana_bond=sqrt((xa(1)-xb(1))**2+(xa(2)-xb(2))**2+(xa(3)-xb(3))**2)
      argos_diana_bond=sqrt(dx(1)**2+dx(2)**2+dx(3)**2)
c
      return
      end
