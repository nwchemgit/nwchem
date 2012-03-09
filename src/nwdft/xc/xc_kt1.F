c
c     this is the gradient correction term by Keal and Tozer,
c     which can be used on its own (as in the KT1 functional),
c     or forms part of the SSB-D functional
c
c     it never includes LDA !!
c     (this is taken care of somewhere else)
c
c     note that even though the energy is assigned here to the
c     exchange, this gradient correction is strictly speaking
c     NOT an exchange functional !
c
c     KT1/KT2 reference:
c     T.W. Keal, D.J. Tozer, JCP 2006, 119, 3015
c     SSB-D reference:
c     M. Swart, M. Sola, F.M. Bickelhaupt, JCP 2009, 131, 094103
c
#ifndef SECOND_DERIV
      Subroutine xc_kt1(tol_rho, fac, rho, delrho, 
     &                     Amat, Cmat, nq, ipol, Ex, qwght,ldew,func)
#else
      Subroutine xc_kt1_d2(tol_rho, fac, rho, delrho, 
     &                     Amat, Amat2, Cmat, Cmat2, nq, ipol,
     &                     Ex, qwght,ldew,func)
#endif
c
C$Id$
c
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
#ifdef SECOND_DERIV
c
c     Second Derivatives of the Exchange Energy Functional
c
      double precision Amat2(nq,NCOL_AMAT2), Cmat2(nq,NCOL_CMAT2)
#endif
c
      double precision DELTA, GAMKT
      Parameter (DELTA = 0.1D0, GAMKT= -0.006d0)
c
c References:
c
c    Keal, Tozer, JCP 119, 3015 (2003), JCP 121, 5654 (2004)
c    Swart, Sola, Bickelhaupt, JCP 131, XXXX (2009)
c    Johnson, Gill & Pople, J. Chem. Phys. 98, 5612 (1993)
c
c***************************************************************************
c
      integer n
      double precision hrho
      double precision rho13, rho43, gamma, g, gdenom, gdenom2
#ifdef SECOND_DERIV
      double precision rho23, rhom23, gdenom3
#endif
c
c     NOTE: the gamma from the KT1 formulation is here called
c           gamkt, gamma is in NWChem reserved for grad**2
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
            if (dsqrt(gamma).gt.tol_rho) then
               gamma = 0.25d0 * gamma
            else
               goto 10
            endif
c
            gdenom = 1.d0 / (rho43 + DELTA)
            gdenom2 = gdenom*gdenom
            g = GAMKT * gamma * gdenom
c
            Ex = Ex + 2.d0*g*qwght(n)*fac
            if (ldew)  func(n) = func(n) + 2.d0*g*fac
            Amat(n,1) = Amat(n,1) - (4d0/3d0)*GAMKT*gamma*rho13*
     &                  fac*gdenom2
            Cmat(n,D1_GAA) = Cmat(n,D1_GAA) + GAMKT*gdenom*fac
c
#ifdef SECOND_DERIV
            rho23 = rho13*rho13
            rhom23 = rho13 / (0.5d0*rho(n,1))
            gdenom3 = gdenom2*gdenom
c
            Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA) + (4d0/3d0)*GAMKT*
     &          gamma*(rho23*7d0/3d0 - DELTA*rhom23/3d0)*gdenom3*fac
            Cmat2(n,D2_RA_GAA) = Cmat2(n,D2_RA_GAA)
     &           - (4d0/3d0)*GAMKT*rho13*gdenom2*fac
c
c      second derivative w.r.t. gamma is zero !
c      (by construction)
c      therefore, nothing added to Cmat2(n,D2_GAA_GAA)
c
#endif
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
            if (dsqrt(gamma).lt.tol_rho) then
               goto 25
            endif
c
            gdenom = 1d0 / (rho43 + DELTA)
            gdenom2 = gdenom*gdenom
            g = GAMKT * gamma * gdenom
c
            Ex = Ex + g*qwght(n)*fac
            if (ldew)  func(n) = func(n) + g*fac
            Amat(n,1) = Amat(n,1) - (4d0/3d0)*GAMKT*gamma*rho13*
     &                 gdenom2*fac
            Cmat(n,D1_GAA) = Cmat(n,D1_GAA) + GAMKT*gdenom*fac
c
#ifdef SECOND_DERIV
            rho23 = rho13*rho13
            rhom23 = rho13 / rho(n,2)
            gdenom3 = gdenom2*gdenom
c
            Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA) + (4d0/3d0)*GAMKT*
     &          gamma*(rho23*7d0/3d0 - DELTA*rhom23/3d0)*gdenom3*fac
            Cmat2(n,D2_RA_GAA) = Cmat2(n,D2_RA_GAA)
     &           - (4d0/3d0)*GAMKT*rho13*gdenom2*fac
c
c      second derivative w.r.t. gamma is zero !
c      (by construction)
c      therefore, nothing added to Cmat2(n,D2_GAA_GAA)
#endif
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
            if (dsqrt(gamma).lt.tol_rho) then
               goto 20
            endif
c           
            gdenom = 1d0 / (rho43 + DELTA)
            gdenom2 = gdenom*gdenom 
            g = GAMKT * gamma * gdenom
c           
            Ex = Ex + g*qwght(n)*fac
            if (ldew)  func(n) = func(n) + g*fac
            Amat(n,2) = Amat(n,2) - (4d0/3d0)*GAMKT*gamma*rho13*
     &            fac*gdenom2
            Cmat(n,D1_GBB) = Cmat(n,D1_GBB) + GAMKT*gdenom*fac
c           
#ifdef SECOND_DERIV
            rho23 = rho13*rho13
            rhom23 = rho13 / rho(n,3)
            gdenom3 = gdenom2*gdenom
c
            Amat2(n,D2_RB_RB) = Amat2(n,D2_RB_RB) + (4d0/3d0)*GAMKT*
     &          gamma*(rho23*7d0/3d0 - DELTA*rhom23/3d0)*gdenom3*fac
            Cmat2(n,D2_RB_GBB) = Cmat2(n,D2_RB_GBB)
     &           - (4d0/3d0)*GAMKT*rho13*gdenom2*fac
c
c      second derivative w.r.t. gamma is zero !
c      (by construction)
c      therefore, nothing added to Cmat2(n,D2_GAA_GAA)
c
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
#include "xc_kt1.F"
#endif
