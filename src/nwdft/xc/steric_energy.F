c     Calculates the steric energy as defined in the reference below
c     for a set of grid points and integration weights

      subroutine steric_energy(tol_rho, fac, rho, 
     &       delrho, nq, qwght, ipol, energy)
c
      implicit none
c      
#include "dft2drv.fh"
c
      double precision tol_rho, fac
      integer nq, ipol
      double precision qwght(nq) ! integration weights
      double precision energy
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
      double precision BETA
      Parameter (BETA = 0.05D0)
c
c     References:
c     S.B. Liu, J. Chem. Phys. 126, 244103 (2007).
c     S.B. Liu, N. Govind, J. Phys. Chem. A. 112, 6690 (2008).
c
c***************************************************************************
c
      integer n
      !double precision arcsinh
      !double precision rho13, rho43, gamma, x, g, gdenom
      double precision gamma, x
c
c     arcsinh(x)=log(x+dsqrt(1d0+x*x))
c
      if (ipol.eq.1) then
c
c        ======> SPIN-RESTRICTED <======
c
         energy = 0.d0 
         do 10 n = 1, nq
            if (rho(n,1).lt.tol_rho) goto 10
c
c           Spin alpha:
c
c           rho13 = (0.5d0*rho(n,1))**(1.d0/3.d0)
c           rho43 = rho13**4 
            gamma = delrho(n,1,1)*delrho(n,1,1) +
     &              delrho(n,2,1)*delrho(n,2,1) +
     &              delrho(n,3,1)*delrho(n,3,1)
            if (gamma.gt.tol_rho)then
               x = gamma/rho(n,1)
            else
               x = 0d0
            endif
            energy = energy + x*qwght(n)
 10      continue
c
      else
c
c        ======> SPIN-UNRESTRICTED <======
c
         energy = 0.d0 
         do 20 n = 1, nq
            if (rho(n,1).lt.tol_rho) goto 20
            if (rho(n,2).lt.tol_rho) goto 25
c
c           Spin alpha:
c
c           rho13 = abs(rho(n,2))**(1.d0/3.d0)*sign(1d0,rho(n,2))
c           rho43 = rho13**4 
            gamma = delrho(n,1,1)*delrho(n,1,1) +
     &              delrho(n,2,1)*delrho(n,2,1) +
     &              delrho(n,3,1)*delrho(n,3,1)
            if (gamma.gt.tol_rho)then
               x = gamma/rho(n,2)
            else
               x = 0d0
            endif
            energy = energy + x*qwght(n)
 25         continue
c
c           Spin beta:
c
            if (rho(n,3).lt.tol_rho) goto 20
c
c           rho13 = abs(rho(n,3))**(1.d0/3.d0)*sign(1d0,rho(n,3))
c           rho43 = rho13**4 
            gamma = delrho(n,1,2)*delrho(n,1,2) +
     &              delrho(n,2,2)*delrho(n,2,2) +
     &              delrho(n,3,2)*delrho(n,3,2)
            if (gamma.gt.tol_rho)then
               x = gamma/rho(n,3)
            else
               x = 0d0
            endif
            energy = energy + x*qwght(n)
 20      continue
c
      endif
c
      energy = energy/8.d0

      return
      end
c $Id$
