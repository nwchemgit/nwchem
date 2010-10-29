c
c $Id$
c
      subroutine hnd_whermt
      implicit none

#include "hnd_whermt.fh"
c
c     From Numerical Recipes: General formulation of abscissas (in h)
c     and weights (in w) of Gauss-Hermite quadrature
c
      double precision eps, pim4
      integer i,ioff,its,j,m,iang,maxit
      parameter (eps=1.d-15,pim4=.7511255444649425d0,maxit=20)
      double precision p1,p2,p3,pp,z,z1
c
      ioff = 0
      do iang =  1, 15
         m = (iang+1)/2
         do i = 1, m
            if (i.eq.1) then
                z = sqrt(float(2*iang+1))-1.85575*(2*iang+1)**(-.16667)
            elseif (i.eq.2) then
                z = z-1.14*iang**.426/z
            elseif (i.eq.3) then
                z = 1.86*z-.86*h(ioff+1)
            elseif (i.eq.4) then
                z = 1.91*z-.91*h(ioff+2)
            else
                z = 2.00*z-h(ioff+i-2)
            endif
            do its=1,maxit
                 p1 = pim4
                 p2 = 0.0d0
                 do j = 1, iang 
                   p3 = p2
                   p2 = p1
                   p1 = z*sqrt(2.0d0/dble(j))*p2
     &                   -sqrt(dble(j-1)/dble(j))*p3
                 enddo
                 pp = sqrt(2.0d0*iang)*p2
                 z1 = z
                 z = z1-p1/pp
                 if (abs(z-z1).le.eps) goto 1
            enddo
            stop 'too many iterations'
  1         h(ioff+i) =  z
            h(ioff+iang+1-i) = -z
            w(ioff+i) = 2.0d0/(pp*pp)
            w(ioff+iang+1-i) = w(ioff+i)
         enddo
         do i = 1, m
            h(ioff+i) = -h(ioff+i)
            h(ioff+iang+1-i) = -h(ioff+iang+1-i)
         enddo
         ioff = ioff + iang
      enddo
      return
      end
