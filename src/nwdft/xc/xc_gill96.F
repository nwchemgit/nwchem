#ifndef SECOND_DERIV
      Subroutine xc_gill96(tol_rho, fac, lfac, nlfac, rho, delrho, 
     &                     Amat, Cmat, nq, ipol, Ex, qwght,ldew,func)
#else
#include "dft2drv.fh"
      Subroutine xc_gill96_d2(tol_rho, fac, lfac, nlfac, rho, delrho, 
     &                        Amat, Amat2, Cmat, Cmat2, nq, ipol, Ex,
     &                        qwght,ldew,func)
#endif
c
C$Id$
c
      implicit none
c      
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
#ifdef SECOND_DERIV
      double precision Amat2(nq,NCOL_AMAT2), Cmat2(nq,NCOL_CMAT2)
#endif
c
      double precision BETA, ONE3, FOUR3
#ifdef SECOND_DERIV
      double precision SEVEN3
#endif
      Parameter (BETA = 1d0/137d0)
      Parameter (ONE3 = 1d0/3d0, FOUR3 = 4d0/3d0)
#ifdef SECOND_DERIV
      Parameter (SEVEN3 = 7d0/3d0)
#endif
c
c References:
c
c    Gill , Mol. Phys. 89, 433 (1996)
c
c***************************************************************************
c
      integer n
      double precision C, rhoval, rho13, rho43, gamma, x, d1x(2),
     &     g, d1g
#ifdef SECOND_DERIV
      double precision d2x(3), d2g
#endif
c
c
c     Uniform electron gas constant
c
      C = -(1.5d0)*(0.75d0/acos(-1d0))**(ONE3)
c
      if (ipol.eq.1) then
c
c        ======> SPIN-RESTRICTED <======
c
         do 10 n = 1, nq
            if (rho(n,1).lt.tol_rho) goto 10
            rhoval = 0.5d0*rho(n,1)
c
c           Spin alpha:
c
            rho13 = rhoval**ONE3
            rho43 = rho13*rhoval
c     Include factor of 4/3 in rho13 since it always appears with it
            rho13 = FOUR3*rho13
c
            if (lfac) then
               Ex = Ex + 2d0*rho43*C*qwght(n)*fac
               if(ldew)func(n) = func(n) + 2.d0*rho43*C*fac
               Amat(n,1) = Amat(n,1) + rho13*C*fac
#ifdef SECOND_DERIV
               Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA)
     &              + ONE3*(rho13/rhoval)*C*fac
#endif
            endif
c
            gamma = delrho(n,1,1)*delrho(n,1,1) +
     &              delrho(n,2,1)*delrho(n,2,1) +
     &              delrho(n,3,1)*delrho(n,3,1)
            if (nlfac.and.dsqrt(gamma).gt.tol_rho)then
               gamma = 0.25d0*gamma
               x = dsqrt(gamma)/rho43
               d1x(1) = -FOUR3*x/rhoval
               d1x(2) = 0.5d0*x/gamma
               g = -BETA*x*sqrt(x)
               d1g = -1.5d0*BETA*sqrt(x)
c
               Ex = Ex + 2d0*rho43*g*qwght(n)*fac
               if(ldew)func(n) = func(n) + 2.d0*rho43*g*fac
               Amat(n,1) = Amat(n,1) + (rho13*g+rho43*d1g*d1x(1))*fac
               Cmat(n,1) = Cmat(n,1) + rho43*d1g*d1x(2)*fac
#ifdef SECOND_DERIV
               d2g = 0.5d0*d1g/x
               d2x(1) = -SEVEN3*d1x(1)/rhoval
               d2x(2) = -FOUR3*d1x(2)/rhoval
               d2x(3) = -0.5d0*d1x(2)/gamma
               Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA)
     &              + ONE3*(rho13/rhoval)*g*fac
     &              + 2.d0*rho13*d1g*d1x(1)*fac
     &              + rho43*(d2g*d1x(1)*d1x(1)+d1g*d2x(1))*fac
               Cmat2(n,D2_RA_GAA) = Cmat2(n,D2_RA_GAA)
     &              + rho13*d1g*d1x(2)*fac
     &              + rho43*(d2g*d1x(1)*d1x(2)+d1g*d2x(2))*fac
               Cmat2(n,D2_GAA_GAA) = Cmat2(n,D2_GAA_GAA)
     &              + rho43*(d2g*d1x(2)*d1x(2)+d1g*d2x(3))*fac
#endif
            endif
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
            rhoval = rho(n,2)
            rho13 = rhoval**ONE3
            rho43 = rho13*rhoval
c     Include factor of 4/3 in rho13 since it always appears with it
            rho13 = FOUR3*rho13
c
            if (lfac) then
               Ex = Ex + rho43*C*qwght(n)*fac
               if(ldew)func(n) = func(n) + rho43*C*fac
               Amat(n,1) = Amat(n,1) + rho13*C*fac
#ifdef SECOND_DERIV
               Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA)
     &              + ONE3*(rho13/rhoval)*C*fac
#endif
            endif
c
            gamma = delrho(n,1,1)*delrho(n,1,1) +
     &              delrho(n,2,1)*delrho(n,2,1) +
     &              delrho(n,3,1)*delrho(n,3,1)
            if (nlfac.and.dsqrt(gamma).gt.tol_rho)then
               x = dsqrt(gamma)/rho43
               d1x(1) = -FOUR3*x/rhoval
               d1x(2) = 0.5d0*x/gamma
               g = -BETA*x*sqrt(x)
               d1g = -1.5d0*BETA*sqrt(x)
c
               Ex = Ex + rho43*g*qwght(n)*fac
               if(ldew)func(n) = func(n) + rho43*g*fac
               Amat(n,1) = Amat(n,1) + (rho13*g+rho43*d1g*d1x(1))*fac
               Cmat(n,1) = Cmat(n,1) + rho43*d1g*d1x(2)*fac
#ifdef SECOND_DERIV
               d2g = 0.5d0*d1g/x
               d2x(1) = -SEVEN3*d1x(1)/rhoval
               d2x(2) = -FOUR3*d1x(2)/rhoval
               d2x(3) = -0.5d0*d1x(2)/gamma
               Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA)
     &              + ONE3*(rho13/rhoval)*g*fac
     &              + 2.d0*rho13*d1g*d1x(1)*fac
     &              + rho43*(d2g*d1x(1)*d1x(1)+d1g*d2x(1))*fac
               Cmat2(n,D2_RA_GAA) = Cmat2(n,D2_RA_GAA)
     &              + rho13*d1g*d1x(2)*fac
     &              + rho43*(d2g*d1x(1)*d1x(2)+d1g*d2x(2))*fac
               Cmat2(n,D2_GAA_GAA) = Cmat2(n,D2_GAA_GAA)
     &              + rho43*(d2g*d1x(2)*d1x(2)+d1g*d2x(3))*fac
#endif
            endif
c
c
 25         continue
c
c           Spin beta:
c
            if (rho(n,3).lt.tol_rho) goto 20
c
            rhoval = rho(n,3)
            rho13 = rhoval**ONE3
            rho43 = rho13*rhoval
c     Include factor of 4/3 in rho13 since it always appears with it
            rho13 = FOUR3*rho13
c
            if (lfac) then
               Ex = Ex + rho43*C*qwght(n)*fac
               if(ldew)func(n) = func(n) + rho43*C*fac
               Amat(n,2) = Amat(n,2) + rho13*C*fac
#ifdef SECOND_DERIV
               Amat2(n,D2_RB_RB) = Amat2(n,D2_RB_RB)
     &              + ONE3*(rho13/rhoval)*C*fac
#endif
            endif
c
            gamma = delrho(n,1,2)*delrho(n,1,2) +
     &              delrho(n,2,2)*delrho(n,2,2) +
     &              delrho(n,3,2)*delrho(n,3,2)
            if (nlfac.and.dsqrt(gamma).gt.tol_rho)then
               x = dsqrt(gamma)/rho43
               d1x(1) = -FOUR3*x/rhoval
               d1x(2) = 0.5d0*x/gamma
               g = -BETA*x*sqrt(x)
               d1g = -1.5d0*BETA*sqrt(x)
c
               Ex = Ex + rho43*g*qwght(n)*fac
               if(ldew)func(n) = func(n) + rho43*g*fac
               Amat(n,2) = Amat(n,2) + (rho13*g+rho43*d1g*d1x(1))*fac
               Cmat(n,3) = Cmat(n,3) + rho43*d1g*d1x(2)*fac
#ifdef SECOND_DERIV
               d2g = 0.5d0*d1g/x
               d2x(1) = -SEVEN3*d1x(1)/rhoval
               d2x(2) = -FOUR3*d1x(2)/rhoval
               d2x(3) = -0.5d0*d1x(2)/gamma
               Amat2(n,D2_RB_RB) = Amat2(n,D2_RB_RB)
     &              + ONE3*(rho13/rhoval)*g*fac
     &              + 2.d0*rho13*d1g*d1x(1)*fac
     &              + rho43*(d2g*d1x(1)*d1x(1)+d1g*d2x(1))*fac
               Cmat2(n,D2_RB_GBB) = Cmat2(n,D2_RB_GBB)
     &              + rho13*d1g*d1x(2)*fac
     &              + rho43*(d2g*d1x(1)*d1x(2)+d1g*d2x(2))*fac
               Cmat2(n,D2_GBB_GBB) = Cmat2(n,D2_GBB_GBB)
     &              + rho43*(d2g*d1x(2)*d1x(2)+d1g*d2x(3))*fac
#endif
            endif
c
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
#include "xc_gill96.F"
#endif
