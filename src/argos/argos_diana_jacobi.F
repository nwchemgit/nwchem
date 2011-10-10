      subroutine argos_diana_jacobi(a,n,na,d,v,nrot)
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
      stop 'argos_diana_jacobi: maximum iterations reached'
c
      return
      end
 
c $Id$
