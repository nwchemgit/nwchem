      subroutine hnd_droot
* $id: hess_quad.f,v 1.5 1998/10/01 20:47:26 windus exp $
      implicit none
#if defined(LINUX) || defined(CRAY) || defined(WIN32)
      double precision ff, r, ww, c, s, a, rt, pt5, zero, one, four
      double precision x, wsum, dum, root, poly
#else
      real*16 ff, r, ww, c, s, a, rt, pt5, zero, one, four
      real*16 x, wsum, dum, root, poly
#endif
#include "hnd_rys.fh"
      integer n, n1, nn, i, j, k, k1, jmax, j1, m
      common/hnd_rysff/ff(19)
      common/hnd_rysrw/r(9,9),ww(9,9)
      dimension c(10,10),s(10,10),a(10),rt(10)
#if defined(LINUX) || defined(CRAY) || defined(WIN32)
      data pt5,zero,one,four /0.5d+00,0.0d+00,1.0d+00,4.0d+00/
#else
      data pt5,zero,one,four /0.5q+00,0.0q+00,1.0q+00,4.0q+00/
#endif
c
c     this version uses christoffel formula for weights.
c     ith root of the jth rys polynomial is returned in r(i,j) with
c     the corresponding weight factor in ww(i,j).   j=1,2,...,n
c
      n=mroots
      x=yy
      if(n.lt.2) n=2
      n1=n+1
      nn=n+n
      call hnd_rysfun(x,nn)
      do 10 i=1,n1
      do 10 j=1,n1
   10 s(i,j)=ff(i+j-1)
      call hnd_ryssmt(c,s,n1)
      do 20 i=1,n
      do 20 j=1,i
      ww(i,j)= zero
   20 r(i,j)= zero
      wsum=ff(1)
      ww(1,1)=wsum
      r(1,1)=ff(2)/wsum
      dum= sqrt(c(2,3)**2-four *c(1,3)*c(3,3))
      r(1,2)=   pt5 *(-c(2,3)-dum)/c(3,3)
      r(2,2)=   pt5 *(-c(2,3)+dum)/c(3,3)
      if(n.eq.2) go to 70
      do 25 i=3,n1
   25 rt(i)=  one
      rt(1)=r(1,2)
      rt(2)=r(2,2)
      do 60 k=3,n
      k1=k+1
      do 30 i=1,k1
   30 a(i)=c(i,k1)
      call hnd_rysnod(a,rt,k)
      do 50 i=1,k
   50 r(i,k)=rt(i)
   60 continue
   70 do 150 k=2,n
      jmax=k-1
      do 150 i=1,k
      root=r(i,k)
      dum=  one  /ff(1)
      do 110 j=1,jmax
      j1=j+1
      poly=c(j1,j1)
      do 100 m=1,j
  100 poly=poly*root+c(j1-m,j1)
  110 dum=dum+poly*poly
  150 ww(i,k)=  one  /dum
      do 160 k=1,mroots
      dum=r(k,mroots)
      u9(k)=dum/(  one  -dum)
  160 w9(k)=ww(k,mroots)
      return
      end
c
      subroutine hnd_rysfun(x,n)
      implicit none
#if defined(LINUX) || defined(CRAY) || defined(WIN32)
      double precision ff, pt5, one, two, x
      double precision tol, e, xx, facmin, term 
      double precision a, sum, fac, t, s, tmax
#else
      real*16 ff, pt5, one, two, x
      real*16 tol, e, xx, facmin, term 
      real*16 a, sum, fac, t, s, tmax
#endif
      integer n, m
      common/hnd_rysff/ff(19)
#if defined(LINUX)  || defined(CRAY) || defined(WIN32)
      data pt5,one,two /0.5d+00,1.0d+00,2.0d+00/
      tol=1.0d-12
      e=0.5409855304296342219319112d-78
      xx=x+x
      facmin=xx
      if(facmin.lt.2*180.2160d+00) e= exp(-x)
      if(facmin.gt.   80.0000d+00) go to 100
#else
      data pt5,one,two /0.5q+00,1.0q+00,2.0q+00/
      tol=1.0q-14
      e=0.5409855304296342219319112q-78
      xx=x+x
      facmin=xx
      if(facmin.lt.2*180.2160q+00) e= exp(-x)
      if(facmin.gt.   80.0000q+00) go to 100
#endif
      term=one
      sum =one
      fac=n
      fac=fac+pt5
   10 fac=fac+one
      term=term*x/fac
      sum=sum+term
      if(fac.le.facmin) go to 10
      t=term
      s=sum
      if(t.gt.s*tol) go to 10
      fac=n+n+1
      ff(n+1)=sum*e/fac
      m=n-1
      fac=m+m+1
   20 if(m.lt.0) return
      ff(m+1)=(e+xx*ff(m+2))/fac
      m=m-1
      fac=fac-two
      go to 20
c
c     use asymptotic expansion for large arguments.
c
#if defined(LINUX) || defined(CRAY) || defined(WIN32)
  100 a= sqrt(.7853981633974483096156608d+00/x)
#else
  100 a= sqrt(.7853981633974483096156608q+00/x)
#endif
      tmax=a*tol/e
      term=one/xx
      sum=term
      fac=one
  110 fac=fac-two
      term=fac*term/xx
      sum=term+sum
      t=term
      if( abs(t).gt.tmax) go to 110
      ff(1)=a-e*sum
      fac=-one
      m=0
  120 if(m.eq.n) return
      m=m+1
      fac=fac+two
      ff(m+1)=(fac*ff(m)-e)/xx
      go to 120
      end
c
      subroutine hnd_ryssmt(c,s,n)
      implicit none
#if defined(LINUX) || defined(CRAY) || defined(WIN32)
      double precision c, s, v, y, zero, one
      double precision fac, dot
#else
      real*16 c, s, v, y, zero, one
      real*16 fac, dot
#endif
      integer n, i, j, k, kmax
      dimension c(10,10),s(10,10),v(10),y(10)
#if defined(LINUX) || defined(CRAY) || defined(WIN32)
      data zero,one /0.0d+00,1.0d+00/
#else
      data zero,one /0.0q+00,1.0q+00/
#endif
c
c     routine returns an n by n triangular matrix c such that
c     c(transpose)sc=i,  where i is an n by n identity matrix.
c
      do 10 i=1,n
      do 10 j=1,i
   10 c(i,j)= zero
      do 100 j=1,n
      kmax=j-1
      fac=s(j,j)
      if(kmax.eq.0) go to 60
      do 20 k=1,kmax
      v(k)= zero
   20 y(k)=s(k,j)
      do 50 k=1,kmax
      dot= zero
      do 30 i=1,k
   30 dot=c(i,k)*y(i)+dot
      do 40 i=1,k
   40 v(i)=v(i)-dot*c(i,k)
   50 fac=fac-dot*dot
   60 fac=one/ sqrt(fac)
      c(j,j)=fac
      if(kmax.eq.0) go to 100
      do 70 k=1,kmax
   70 c(k,j)=fac*v(k)
  100 continue
      return
      end
c
      subroutine hnd_rysnod(a,rt,k)
      implicit none
#if defined(LINUX) || defined(CRAY) || defined(WIN32)
      double precision a, rt, zero, tol, r1, p1, r2, p2
      double precision prod, r3, p3, r4, p4, r5, p5, r6, p6
      double precision r, dr, delta
#else
      real*16 a, rt, zero, tol, r1, p1, r2, p2
      real*16 prod, r3, p3, r4, p4, r5, p5, r6, p6
      real*16 r, dr, delta
#endif
#include "stdio.fh"
      integer k, k1, m, i
      dimension a(10),rt(10)
#if defined(LINUX) || defined(CRAY) || defined(WIN32)
      data zero /0.0d+00/
#else
      data zero /0.0q+00/
#endif
c
c     routine returns rt(i) the ith root of a polynomial of order
c     k whose mth coefficient is stored in a(m+1).  it is assumed that
c     the initial values in rt bracket the final values.
c
#if defined(LINUX)  || defined(CRAY) || defined(WIN32)
      tol=1.0d-09
#else
      tol=1.0q-11
#endif
      k1=k+1
      r2= zero
      p2=a(1)
      do 100 m=1,k
      r1=r2
      p1=p2
      r2=rt(m)
      p2=a(k1)
      do 10 i=1,k
   10 p2=p2*r2+a(k1-i)
      prod=p1*p2
      if(prod.lt. zero) go to 20
      write(luout,15) m,k
   15 format(//,' root number ',i3,' was not found for polynomial',
     1 ' of order ',i3,//)
      call errquit('hnd_rysnod: root not found', 555)
   20 r5=r1
      p5=p1
      r6=r2
      p6=p2
   30 r3=r5
      p3=p5
      r4=r6
      p4=p6
      r =(r3*p4-r4*p3)/(p4-p3)
      dr=r4-r3
      delta=dr
      if( abs(delta).lt.tol) go to 90
#if defined(LINUX) || defined(CRAY) || defined(WIN32)
      dr=0.0625d+00*dr
#else
      dr=0.0625q+00*dr
#endif
      r5=r-dr
      if(r5.lt.r3) r5=r3
      r6=r+dr
      if(r6.gt.r4) r6=r4
      p5=a(k1)
      p6=p5
      do 40 i=1,k
      p5=p5*r5+a(k1-i)
   40 p6=p6*r6+a(k1-i)
   45 prod=p5*p6
      if(prod.lt. zero) go to 30
      prod=p3*p5
      if(prod.gt. zero) go to 60
#if defined(LINUX) || defined(CRAY) || defined(WIN32)
      r5=0.25d+00*r3+0.75d+00*r5
#else
      r5=0.25q+00*r3+0.75q+00*r5
#endif
      p5=a(k1)
      do 50 i=1,k
   50 p5=p5*r5+a(k1-i)
      go to 45
#if defined(LINUX) || defined(CRAY) || defined(WIN32)
   60 r6=0.25d+00*r4+0.75d+00*r6
#else
   60 r6=0.25q+00*r4+0.75q+00*r6
#endif
      p6=a(k1)
      do 70 i=1,k
   70 p6=p6*r6+a(k1-i)
      go to 45
   90 rt(m)=r
  100 continue
      return
      end
