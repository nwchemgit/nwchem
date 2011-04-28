c
C$Id$
c
      subroutine xc_optx(tol_rho, fac, p,  rho, delrho, 
     &                      Amat, Cmat, nq, ipol, Ex, qwght,ldew,func)
      implicit none
c      
#include "dft2drv.fh"
c
      double precision tol_rho, fac, Ex
      integer nq, ipol
      logical ldew
      double precision func(*)  ! value of the functional [output]
c
c     Charge Density
c
      double precision rho(nq,ipol*(ipol+1)/2)
c
c     Charge Density Gradient
c
      double precision delrho(nq,3,ipol)
c
c     Quadrature Weights
c
      double precision qwght(nq)
c
c     Sampling Matrices for the XC Potential
c
      double precision Amat(nq,ipol), Cmat(nq,*)
c
      integer p ! [in]
c
c
c References:
c
c    Becke,  (1986)
c    Handy NC, Cohen AJ, Mol Phys 99 (5): 403-412 MAR 2001
c    idem, Mol Phys 99 (7); 607-615 2001
c
c***************************************************************************
c
      integer n
      double precision rho13, rho43, gamma, x, g,  dg,
     &     t, hrho
      double precision gamma86
      Parameter (gamma86=0.006d0)
      double precision ux,uxp,gx
      ux(x,gx)=gx*x*x/(1d0+gx*x*x)
      uxp(x,gx)=gx*x*2d0/(1d0+gx*x*x)**2
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
               x = sqrt(gamma) / rho43
            else
               x = 0d0
            endif
c
            g = -ux(x,gamma86)**p
            dg = -p*ux(x,gamma86)**(p-1)*uxp(x,gamma86)
c
c
            Ex = Ex + 2d0*rho43*g*qwght(n)*fac
            if(ldew)func(n) = func(n) + 2.d0*rho43*g*fac
            Amat(n,1) = Amat(n,1) + (4d0/3d0)*rho13*(g-x*dg)*fac
c
            if (x.gt.tol_rho) then
               t = 0.5d0 * dg / sqrt(gamma) * fac
               Cmat(n,D1_GAA) = Cmat(n,D1_GAA) + t
            endif
c
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
            g = -ux(x,gamma86)**p
            dg = -p*ux(x,gamma86)**(p-1)*uxp(x,gamma86)
c
c
            Ex = Ex + rho43*g*qwght(n)*fac
            if (ldew)func(n) = func(n) + rho43*g*fac
            Amat(n,1) = Amat(n,1) + (4d0/3d0)*rho13*(g-x*dg)*fac
c     
            if (x.gt.tol_rho) then
               t = dg / sqrt(gamma) * fac
               Cmat(n,D1_GAA) = Cmat(n,D1_GAA) + t * 0.5d0
            endif
c
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
            g = -ux(x,gamma86)**p
            dg = -p*ux(x,gamma86)**(p-1)*uxp(x,gamma86)
c
c

            Ex = Ex + rho43*g*qwght(n)*fac
            if (ldew)func(n) = func(n) +rho43*g*fac
            Amat(n,2) = Amat(n,2) + (4d0/3d0)*rho13*(g-x*dg)*fac
c     
            if (x.gt.tol_rho) then
               t = dg / sqrt(gamma) * fac
               Cmat(n,D1_GBB) = Cmat(n,D1_GBB) + t * 0.5d0
            endif
c
c
 20      continue
c
      endif
c
      return
      end
