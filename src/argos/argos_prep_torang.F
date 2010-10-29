      logical function argos_prep_torang(x,mx,ia,ja,ka,la,angle)
c
c $Id$
c
      implicit none
c
      integer mx,ia,ja,ka,la
      real*8 x(3,mx),angle
c
      integer i,j,k
      real*8 xij(3),xkj(3),xkl(3),xm(3),xn(3)
c
      do 1 i=1,3
      xij(i)=x(i,ia)-x(i,ja)
      xkj(i)=x(i,ka)-x(i,ja)
      xkl(i)=x(i,ka)-x(i,la)
    1 continue
c
      do 2 i=1,3
      j=i+1
      if(j.eq.4) j=1
      k=j+1
      if(k.eq.4) k=1
      xm(i)=xij(j)*xkj(k)-xkj(j)*xij(k)
      xn(i)=xkj(j)*xkl(k)-xkl(j)*xkj(k)
    2 continue
c
      angle=acos((xm(1)*xn(1)+xm(2)*xn(2)+xm(3)*xn(3))/
     + sqrt((xm(1)*xm(1)+xm(2)*xm(2)+xm(3)*xm(3))*
     + (xn(1)*xn(1)+xn(2)*xn(2)+xn(3)*xn(3))))
c
      if(xkj(1)*(xm(2)*xn(3)-xm(3)*xn(2))+
     + xkj(2)*(xm(3)*xn(1)-xm(1)*xn(3))+
     + xkj(3)*(xm(1)*xn(2)-xm(2)*xn(1)).lt.0.0d0) angle=-angle
c
      argos_prep_torang=.true.
      return
      end
