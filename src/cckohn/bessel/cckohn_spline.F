      subroutine cckohn_spline(x,y,n,yp1,ypn,scr,y2)
c
      implicit none
      integer n
      double precision x(n),y(n),yp1,ypn,scr(n),y2(n),p,qn,un,sig
      integer i,k
c
      y2(1)=-.5
      scr(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      do 11 i=2,n-1
         sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
         p=sig*y2(i-1)+2.
         y2(i)=(sig-1.)/p
         scr(i)=(6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))
     &        /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*scr(i-1))/p
 11   continue
      qn=.5
      un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      y2(n)=(un-qn*scr(n-1))/(qn*y2(n-1)+1.)
      do 12 k=n-1,1,-1
         y2(k)=y2(k)*y2(k+1)+scr(k)
 12   continue
      return
      end
c $Id$
