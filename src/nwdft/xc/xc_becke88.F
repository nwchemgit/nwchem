#if defined(FUJITSU_VPP)
!ocl scalar
#endif
#ifndef SECOND_DERIV
      Subroutine xc_becke88(tol_rho, fac, lfac, nlfac, rho, delrho, 
     &                      Amat, Cmat, nq, ipol, Ex, qwght,ldew,func)
#else
      Subroutine xc_becke88_d2(tol_rho, fac, lfac, nlfac, rho, delrho, 
     &                         Amat, Amat2, Cmat, Cmat2, nq, ipol, Ex,
     &                         qwght,ldew,func)
#endif
c
C$Id: xc_becke88.F,v 1.22 2003-03-06 20:52:58 edo Exp $
c
      implicit none
c      
#include "dft2drv.fh"
c
      double precision tol_rho, fac, Ex
      integer nq, ipol
      logical lfac, nlfac,ldew
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
#ifdef SECOND_DERIV
c
c     Second Derivatives of the Exchange Energy Functional
c
      double precision Amat2(nq,NCOL_AMAT2), Cmat2(nq,NCOL_CMAT2)
#endif
c
      double precision BETA
      Parameter (BETA = 0.0042D0)
c
c References:
c
c    Becke, Phys. Rev. A 38, 3098 (1988)
c    Johnson, Gill & Pople, J. Chem. Phys. 98, 5612 (1993)
c
c***************************************************************************
c
      integer n
      double precision arcsinh, darcsinh
      double precision C, rho13, rho43, gamma, x, g, gdenom, dg,
     &     dgdenom, t
#ifdef SECOND_DERIV
      double precision rhom23, d2g, d2gdenom
#endif
c
      arcsinh(x)=log(x+dsqrt(1d0+x*x))
      darcsinh(x)=1d0/dsqrt(1d0+x*x)
c
c     Uniform electron gas constant
c
      C = -(1.5d0)*(0.75d0/acos(-1d0))**(1d0/3d0)
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
            gdenom = 1d0 + 6d0*BETA*x*arcsinh(x)
            g = -BETA*x*x / gdenom
            dgdenom = 6d0*BETA*(arcsinh(x) + x*darcsinh(x))
            dg = BETA*x*(x*dgdenom - 2d0*gdenom) / gdenom**2
c
            if (lfac) then
               Ex = Ex + 2d0*rho43*C*qwght(n)*fac
               if(ldew)func(n) = func(n) + 2.d0*rho43*C*fac
               Amat(n,1) = Amat(n,1) + (4d0/3d0)*rho13*C*fac
            endif
c
            if (nlfac) then
               Ex = Ex + 2d0*rho43*g*qwght(n)*fac
               if(ldew)func(n) = func(n) + 2.d0*rho43*g*fac
               Amat(n,1) = Amat(n,1) + (4d0/3d0)*rho13*(g-x*dg)*fac
            endif
c
            if (x.gt.tol_rho) then
               t = 0.5d0 * dg / sqrt(gamma) * fac
               Cmat(n,D1_GAA) = Cmat(n,D1_GAA) + t
c               Cmat(n,1,1) = Cmat(n,1,1) + t*delrho(n,1,1)
c               Cmat(n,2,1) = Cmat(n,2,1) + t*delrho(n,2,1)
c               Cmat(n,3,1) = Cmat(n,3,1) + t*delrho(n,3,1)
            endif
c
#ifdef SECOND_DERIV
            if(lfac) g = g + C           ! Add local contribution back to g
            rhom23 = rho13 / (0.5d0*rho(n,1))
            d2gdenom = 6d0*BETA*darcsinh(x)*(2d0 - x*x/(x*x+1d0))
            d2g = -2d0*BETA/gdenom + 4d0*BETA*x*dgdenom/gdenom**2
     &           + BETA*x*x*d2gdenom/gdenom**2
     &           - 2d0*BETA*x*x*(dgdenom)**2/gdenom**3
c
            Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA)
     &           + (4d0/9d0)*rhom23*(g-x*dg+4d0*x*x*d2g)*fac
            Cmat2(n,D2_RA_GAA) = Cmat2(n,D2_RA_GAA)
     &           - (4d0/3d0)*(rhom23**2/rho(n,1))*d2g*fac
            if (x.gt.tol_rho) then
               Cmat2(n,D2_GAA_GAA) = Cmat2(n,D2_GAA_GAA)
     &              - 0.25d0*gamma**(-1.5d0)*(dg-x*d2g)*fac
            endif
#endif
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
            gdenom = 1d0 + 6d0*BETA*x*arcsinh(x)
            g = -BETA*x*x / gdenom
            dgdenom = 6d0*BETA*(arcsinh(x) + x*darcsinh(x))
            dg = BETA*x*(x*dgdenom - 2d0*gdenom) / gdenom**2
c
            if (lfac) then
               Ex = Ex + rho43*C*qwght(n)*fac
               if (ldew)func(n) = func(n) + rho43*C*fac
               Amat(n,1) = Amat(n,1) + (4d0/3d0)*rho13*C*fac
            endif
c
            if (nlfac) then
               Ex = Ex + rho43*g*qwght(n)*fac
               if (ldew)func(n) = func(n) + rho43*g*fac
               Amat(n,1) = Amat(n,1) + (4d0/3d0)*rho13*(g-x*dg)*fac
            endif
c
            if (x.gt.tol_rho) then
               t = dg / sqrt(gamma) * fac
               Cmat(n,D1_GAA) = Cmat(n,D1_GAA) + t * 0.5d0
c               Cmat(n,1,1) = Cmat(n,1,1) + t*delrho(n,1,1)
c               Cmat(n,2,1) = Cmat(n,2,1) + t*delrho(n,2,1)
c               Cmat(n,3,1) = Cmat(n,3,1) + t*delrho(n,3,1)
            endif
c
#ifdef SECOND_DERIV
            if (lfac) g = g + C           ! Add local contribution back to g
            rhom23 = rho13 / rho(n,2)
            d2gdenom = 6d0*BETA*darcsinh(x)*(2d0 - x*x/(x*x+1d0))
            d2g = -2d0*BETA/gdenom + 4d0*BETA*x*dgdenom/gdenom**2
     &           + BETA*x*x*d2gdenom/gdenom**2
     &           - 2d0*BETA*x*x*(dgdenom)**2/gdenom**3
c
            Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA)
     &           + (4d0/9d0)*rhom23*(g-x*dg+4d0*x*x*d2g)*fac
            Cmat2(n,D2_RA_GAA) = Cmat2(n,D2_RA_GAA)
     &           - (2d0/3d0)*(rhom23**2/rho(n,2))*d2g*fac
            if (x.gt.tol_rho) then
               Cmat2(n,D2_GAA_GAA) = Cmat2(n,D2_GAA_GAA)
     &              - 0.25d0*gamma**(-1.5d0)*(dg-x*d2g)*fac
            endif
#endif
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
            gdenom = 1d0 + 6d0*BETA*x*arcsinh(x)
            g = -BETA*x*x / gdenom
            dgdenom = 6d0*BETA*(arcsinh(x) + x*darcsinh(x))
            dg = BETA*x*(x*dgdenom - 2d0*gdenom) / gdenom**2
c
            if (lfac) then
               Ex = Ex + rho43*C*qwght(n)*fac
               if (ldew)func(n) = func(n) + rho43*C*fac
               Amat(n,2) = Amat(n,2) + (4d0/3d0)*rho13*C*fac
            endif
c
            if (nlfac) then
               Ex = Ex + rho43*g*qwght(n)*fac
               if (ldew)func(n) = func(n) +rho43*g*fac
               Amat(n,2) = Amat(n,2) + (4d0/3d0)*rho13*(g-x*dg)*fac
            endif
c
            if (x.gt.tol_rho) then
               t = dg / sqrt(gamma) * fac
               Cmat(n,D1_GBB) = Cmat(n,D1_GBB) + t * 0.5d0
c               Cmat(n,1,2) = Cmat(n,1,2) + t*delrho(n,1,2)
c               Cmat(n,2,2) = Cmat(n,2,2) + t*delrho(n,2,2)
c               Cmat(n,3,2) = Cmat(n,3,2) + t*delrho(n,3,2)
            endif
c
#ifdef SECOND_DERIV
            if(lfac) g = g + C           ! Add local contribution back to g
            rhom23 = rho13 / rho(n,3)
            d2gdenom = 6d0*BETA*darcsinh(x)*(2d0 - x*x/(x*x+1d0))
            d2g = -2d0*BETA/gdenom + 4d0*BETA*x*dgdenom/gdenom**2
     &           + BETA*x*x*d2gdenom/gdenom**2
     &           - 2d0*BETA*x*x*(dgdenom)**2/gdenom**3
c
            Amat2(n,D2_RB_RB) = Amat2(n,D2_RB_RB)
     &           + (4d0/9d0)*rhom23*(g-x*dg+4d0*x*x*d2g)*fac
            Cmat2(n,D2_RB_GBB) = Cmat2(n,D2_RB_GBB)
     &           - (2d0/3d0)*(rhom23**2/rho(n,3))*d2g*fac
            if (x.gt.tol_rho) then
               Cmat2(n,D2_GBB_GBB) = Cmat2(n,D2_GBB_GBB)
     &              - 0.25d0*gamma**(-1.5d0)*(dg-x*d2g)*fac
            endif
#endif
c
 20      continue
c
      endif
c
      return
      end
#ifndef SECOND_DERIV
#define SECOND_DERIV
c
c     Compile source again for the 2nd derivative case
c
#include "xc_becke88.F"
#endif
