c
c $Id$
c
      subroutine hnd_wleg
      implicit none

#include "hnd_wleg.fh"
c
c     From Numerical Recipes: General formulation of abscissas (in h)
c     and weights (in w) of Gauss-Legrende quadrature
c
      double precision eps
      parameter (eps=1.d-15)
      integer m,n,i,j
      double precision pi,z,z1,p1,p2,p3,pp
c
      n =gleg_points
      pi=acos(-1.0d0)
      m =(n+1)/2
      do i=1,m
         z=cos(pi*(float(i)-0.25d0)/(float(n)+0.5d0))
 1       p1=1.0d0
         p2=0.0d0
         do j=1,n
            p3=p2
            p2=p1
            p1=((2.0d0*dble(j)-1.0d0)*z*p2-(dble(j)-1.0d0)*p3)/(dble(j))
         enddo
         pp=dble(n)*(z*p1-p2)/(z*z-1.0d0)
         z1=z
         z=z1-p1/pp
         if(abs(z-z1).gt.eps) goto 1
         hl(i)=-z
         hl(n+1-i)=z
         wl(i)=2.0d0/((1.0d0-z*z)*pp*pp)
         wl(n+1-i)=wl(i)
      enddo
      return
      end
