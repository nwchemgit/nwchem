#if defined(FUJITSU_VPP)
!ocl scalar
#endif
      Subroutine xc_lb94(tol_rho, fac, rho, delrho, 
     &                      Amat, nq, ipol)
c
C$Id: xc_lb94.F,v 1.5 2003-04-10 17:49:33 sohirata Exp $
c
      implicit none
c      
#include "dft2drv.fh"
c
      double precision tol_rho, fac
      integer nq, ipol
c
c     Charge Density
c
      double precision rho(nq,ipol*(ipol+1)/2)
c
c     Charge Density Gradient
c
      double precision delrho(nq,3,ipol)
c
c     Sampling Matrices for the XC Potential
c
      double precision Amat(nq,ipol)
c
c
      double precision BETA
      Parameter (BETA = 0.05D0)
c
c References:
c
c    R. van Leeuwen & E. J. Baerends, Phys. Rev. A 49, 2421 (1994).
c
c***************************************************************************
c
      integer n
      double precision arcsinh
      double precision rho13, rho43, gamma, x, g, gdenom
c
      arcsinh(x)=log(x+dsqrt(1d0+x*x))
c
      if (ipol.eq.1) then
c
c        ======> SPIN-RESTRICTED <======
c
         do 10 n = 1, nq
            if (rho(n,1).lt.tol_rho) goto 10
c
c           Spin alpha:
c
            rho13 = (0.5d0*rho(n,1))**(1.d0/3.d0)
            rho43 = rho13**4 
            gamma = delrho(n,1,1)*delrho(n,1,1) +
     &              delrho(n,2,1)*delrho(n,2,1) +
     &              delrho(n,3,1)*delrho(n,3,1)
            if (dsqrt(gamma).gt.tol_rho)then
               gamma = 0.25d0 * gamma
               x = sqrt(gamma) / rho43
            else
               x = 0d0
            endif
c
            gdenom = 1d0 + 3d0*BETA*x*arcsinh(x)
            g = -BETA*x*x / gdenom
c
            Amat(n,1) = Amat(n,1) + rho13*g*(1d0-fac)
c
 10      continue
c
      else
c
c        ======> SPIN-UNRESTRICTED <======
c
         do 20 n = 1, nq
            if (dabs(rho(n,1)).lt.tol_rho) goto 20
            if (dabs(rho(n,2)).lt.tol_rho) goto 25
c
c           Spin alpha:
c
            rho13 = abs(rho(n,2))**(1.d0/3.d0)*sign(1d0,rho(n,2))
            rho43 = rho13**4 
            gamma = delrho(n,1,1)*delrho(n,1,1) +
     &              delrho(n,2,1)*delrho(n,2,1) +
     &              delrho(n,3,1)*delrho(n,3,1)
            if (dsqrt(gamma).gt.tol_rho)then
               x = sqrt(gamma) / rho43
            else
               x = 0d0
            endif
c
            gdenom = 1d0 + 3d0*BETA*x*arcsinh(x)
            g = -BETA*x*x / gdenom
c
            Amat(n,1) = Amat(n,1) + rho13*g*(1d0-fac)
c
 25         continue
c
c           Spin beta:
c
            if (dabs(rho(n,3)).lt.tol_rho) goto 20
c
            rho13 = abs(rho(n,3))**(1.d0/3.d0)*sign(1d0,rho(n,3))
            rho43 = rho13**4 
            gamma = delrho(n,1,2)*delrho(n,1,2) +
     &              delrho(n,2,2)*delrho(n,2,2) +
     &              delrho(n,3,2)*delrho(n,3,2)
            if (dsqrt(gamma).gt.tol_rho)then
               x = sqrt(gamma) / rho43
            else
               x = 0d0
            endif
c
            gdenom = 1d0 + 3d0*BETA*x*arcsinh(x)
            g = -BETA*x*x / gdenom
c
            Amat(n,2) = Amat(n,2) + rho13*g*(1d0-fac)
c
 20      continue
c
      endif
c
      return
      end
