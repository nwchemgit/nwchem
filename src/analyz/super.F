      subroutine super(x,ix,nx,mx,y,w,sdev,ny,my,mod,rms0,rms1)
c
c $Id$
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
      integer nx,mx,ny,my
      integer ix(mx)
      real*8 x(mx,3),y(my,3),w(my),sdev(my),rms0,rms1
      logical mod
c
      integer i,j,k,nr
      real*8 u(3,3),c(4,4),q(4),b(4),v(4,4),wnorm
      real*8 xt(3),yt(3),sd
c
      wnorm=zero
      do 1 j=1,3
      xt(j)=zero
      yt(j)=zero
    1 continue
      do 2 k=1,nx
      do 3 j=1,3
      xt(j)=xt(j)+w(ix(k))*x(k,j)
    3 continue
    2 continue
      call ga_dgop(mag_d08,xt,3,'+')
      do 4 k=1,ny
      wnorm=wnorm+w(k)
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
      rms0=rms0+w(ix(k))*(x(k,j)-xt(j)-y(ix(k),j)+yt(j))**2
    7 continue
    6 continue
      call ga_dgop(mag_d09,rms0,1,'+')
      rms0=sqrt(rms0/wnorm)
c
      do 8 i=1,3
      do 9 j=1,3
      u(i,j)=zero
      do 10 k=1,nx
      u(i,j)=u(i,j)+w(ix(k))*(x(k,i)-xt(i))*(y(ix(k),j)-yt(j))
   10 continue
    9 continue
    8 continue
      call ga_dgop(mag_d10,u,9,'+')
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
      call jacobi(c,4,4,b,v,nr)
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
      sd=w(ix(k))*sd
      sdev(ix(k))=sdev(ix(k))+sd
      rms1=rms1+sd
   14 continue
      call ga_dgop(mag_d11,rms1,1,'+')
      rms1=sqrt(rms1/wnorm)
c
      return
      end
c
      subroutine jacobi(a,n,na,d,v,nrot)
c
c     compute eigenvectors and eigenvalues for real symmetric 
c     matrix using the Jacobi diagonalization
c
      implicit none
c
      integer nmax,nrmax
      parameter(nmax=100)
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
      stop 'jacobi: maximum iterations reached'
c
      return
      end

