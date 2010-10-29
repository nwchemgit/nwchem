      subroutine esp_jacobi(a)
c
c $Id$
c
      implicit none
c
#include "esp_params.fh"
#include "esp_common.fh"
c
      integer nmax
      parameter (nmax=100)
c
      integer n,nrot
      real*8 a(3,3),d(3),v(3,3),b(nmax),z(nmax)
c
      integer i,j,k,l
      real*8 c,sum,tresh,g,h,s,t,tau,theta
c
      n=3
c
      do 1 i=1,n
      do 2 j=1,n
      v(i,j)=zero
    2 continue
      v(i,i)=one
    1 continue
c
      do 3 i=1,n
      b(i)=a(i,i)
      d(i)=b(i)
      z(i)=zero
    3 continue
      nrot=0
      do 4 k=1,nmax
      sum=zero
      do 5 i=1,n-1
      do 6 j=i+1,n
      sum=sum+abs(a(i,j))
    6 continue
    5 continue
      if(sum.eq.0.)return
      if(k.lt.4)then
      tresh=0.2d0*sum/n**2
      else
      tresh=zero
      endif
      do 7 i=1,n-1
      do 8 j=i+1,n
      g=hundrd*abs(a(i,j))
      if((i.gt.4).and.(abs(d(i))+g.eq.abs(d(i)))
     + .and.(abs(d(j))+g.eq.abs(d(j))))then
      a(i,j)=zero
      else if(abs(a(i,j)).gt.tresh)then
      h=d(j)-d(i)
      if(abs(h)+g.eq.abs(h))then
      t=a(i,j)/h
      else
      theta=half*h/a(i,j)
      t=one/(abs(theta)+sqrt(one+theta**2))
      if(theta.lt.zero)t=-t
      endif
      c=1./sqrt(1+t**2)
      s=t*c
      tau=s/(1.+c)
      h=t*a(i,j)
      z(i)=z(i)-h
      z(j)=z(j)+h
      d(i)=d(i)-h
      d(j)=d(j)+h
      a(i,j)=zero
      do 9 l=1,i-1
      g=a(l,i)
      h=a(l,j)
      a(l,i)=g-s*(h+g*tau)
      a(l,j)=h+s*(g-h*tau)
    9 continue
      do 10 l=i+1,j-1
      g=a(i,l)
      h=a(l,j)
      a(i,l)=g-s*(h+g*tau)
      a(l,j)=h+s*(g-h*tau)
   10 continue
      do 11 l=j+1,n
      g=a(i,l)
      h=a(j,l)
      a(i,l)=g-s*(h+g*tau)
      a(j,l)=h+s*(g-h*tau)
   11 continue
      do 20 l=1,n
      g=v(l,i)
      h=v(l,j)
      v(l,i)=g-s*(h+g*tau)
      v(l,j)=h+s*(g-h*tau)
   20 continue
      nrot=nrot+1
      endif
    8 continue
    7 continue
      do 21 i=1,n
      b(i)=b(i)+z(i)
      d(i)=b(i)
      z(i)=zero
   21 continue
    4 continue
c
      call md_abort('esp_jacobi: maximum iterations reached',me)
c
      return
      end

