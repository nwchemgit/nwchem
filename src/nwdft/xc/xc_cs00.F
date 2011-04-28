#if defined(FUJITSU_VPP)
!ocl scalar
#endif
      Subroutine xc_cs00(tol_rho, fac, rho, delrho, 
     &                      Amat, nq, ipol, delta, e_homo)
c
C$Id$
c
      implicit none
c      
#include "dft2drv.fh"
c
      double precision tol_rho, fac, delta, shift
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
c    M. E. Casida & D. R. Salahub, J. Chem. Phys. 113, 8918 (2000).
c
c***************************************************************************
c
      integer n
      double precision arcsinh, p1, p2
      double precision rho13, rho43, gamma, x, g, gdenom, hrho
      double precision e_homo
      arcsinh(x)=log(x+dsqrt(1d0+x*x))
c
c     Zhan-Nichols-Dixon shift
c     CAUTION: this shift works only with B3LYP/aug-cc-pVTZ+
c
      if (delta.gt.1.0d90) then
         shift = -0.2332D0 * e_homo + 0.0116D0
      else
         shift = delta
      endif
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
            hrho  = 0.5d0*rho(n,1)
            rho13 = hrho**(1.d0/3.d0)
            rho43 = rho13*hrho
            gamma = delrho(n,1,1)*delrho(n,1,1) +
     &              delrho(n,2,1)*delrho(n,2,1) +
     &              delrho(n,3,1)*delrho(n,3,1)
            if (dsqrt(gamma).gt.tol_rho)then
               gamma = 0.25d0 * gamma
               x = dsqrt(gamma) / rho43
            else
               x = 0d0
            endif
c
            gdenom = 1d0 + 3d0*BETA*x*arcsinh(x)
            g = -BETA*x*x / gdenom
c
            p1 = Amat(n,1) + rho13*g*(1d0-fac)
            p2 = Amat(n,1) - shift
            Amat(n,1) = max(p1, p2)
c
 10      continue
c
      else
c
c        ======> SPIN-UNRESTRICTED <======
c
         do 20 n = 1, nq
            if (rho(n,1).lt.tol_rho) goto 20
            if (rho(n,2).lt.tol_rho) goto 25
c
c           Spin alpha:
c
            rho13 = rho(n,2)**(1.d0/3.d0)
            rho43 = rho13*rho(n,2)
            gamma = delrho(n,1,1)*delrho(n,1,1) +
     &              delrho(n,2,1)*delrho(n,2,1) +
     &              delrho(n,3,1)*delrho(n,3,1)
            if (dsqrt(gamma).gt.tol_rho)then
               x = dsqrt(gamma) / rho43
            else
               x = 0d0
            endif
c
            gdenom = 1d0 + 3d0*BETA*x*arcsinh(x)
            g = -BETA*x*x / gdenom
c
            p1 = Amat(n,1) + rho13*g*(1d0-fac)
            p2 = Amat(n,1) - shift
            Amat(n,1) = max(p1, p2)
c
 25         continue
c
c           Spin beta:
c
            if (rho(n,3).lt.tol_rho) goto 20
c
            rho13 = rho(n,3)**(1.d0/3.d0)
            rho43 = rho13*rho(n,3)
            gamma = delrho(n,1,2)*delrho(n,1,2) +
     &              delrho(n,2,2)*delrho(n,2,2) +
     &              delrho(n,3,2)*delrho(n,3,2)
            if (dsqrt(gamma).gt.tol_rho)then
               x = dsqrt(gamma) / rho43
            else
               x = 0d0
            endif
c
            gdenom = 1d0 + 3d0*BETA*x*arcsinh(x)
            g = -BETA*x*x / gdenom
c
            p1 = Amat(n,2) + rho13*g*(1d0-fac)
            p2 = Amat(n,2) - shift
            Amat(n,2) = max(p1, p2)
c
 20      continue
c
      endif
c
      return
      end
