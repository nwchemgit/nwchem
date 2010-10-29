      subroutine writel(p,newbas,oprint)
C$Id$
      implicit none
c
      logical oprint
      real *8 p(*)
      integer newbas
c
      integer ind, i, j
c
      integer m5, m, n
c
      ind(i,j)=max0(i,j)*(max0(i,j)-1)/2+min0(i,j)
c
      m5=12
      if(oprint) m5=8
      m=1
      n=m5
   6  if(newbas.lt.m)return
      if(n.gt.newbas)n=newbas
      if(.not.oprint)write(6,200)(i,i=m,n)
      if(oprint)write(6,100)(i,i=m,n)
100   format(//3x,8i14)
200   format(//12i9)
      write(6,101)
101   format(/)
      do 1 j=1,newbas
      if(oprint)write(6,102)j,(p(ind(i,j)),i=m,n)
      if(.not.oprint)write(6,202)j,(p(ind(i,j)),i=m,n)
 1    continue
102   format(7x,i3,8f14.7)
202   format(1x,i3,12f9.4)
      m=m+m5
      n=n+m5
      goto 6
      end
