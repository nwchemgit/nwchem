      subroutine gauleg(x1,x2,x,w,n)
*
* $Id$
*
      implicit double precision  (a-h,o-z)
      double precision  x1,x2,x(n),w(n)
      parameter (eps=3.d-14)
      m=(n+1)/2
      xm=0.5d0*(x2+x1)
      xl=0.5d0*(x2-x1)
      pi = 4.0d0*atan(1.0d0)
      do 12 i=1,m
        z=cos(pi*(i-.25d0)/(n+.5d0))
1       continue
          p1=1.d0
          p2=0.d0
          do 11 j=1,n
            p3=p2
            p2=p1
            p1=((2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/j
11        continue
          pp=n*(z*p1-p2)/(z*z-1.d0)
          z1=z
          z=z1-p1/pp
        if(abs(z-z1).gt.eps)go to 1
        x(i)=xm-xl*z
        x(n+1-i)=xm+xl*z
        w(i)=2.d0*xl/((1.d0-z*z)*pp*pp)
        w(n+1-i)=w(i)
12    continue
      return
      end
