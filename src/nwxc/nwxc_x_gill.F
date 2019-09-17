#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#if !defined(NWAD_PRINT)
C> \ingroup nwxc
C> @{
C>
C> \file nwxc_x_gill.F
C> The gradient-corrected exchange functional of Gill
C>
C> @}
#endif
#endif
C> \ingroup nwxc_priv
C>
C> \brief Evaluates the Gill exchange functional
C>
C> Evaluates the Gill exchange functional which is a simplified
C> functional that yields results similar to Becke'88 [1].
C>
C> ### References ###
C>
C> [1] P. M. W. Gill,
C>     "A new gradient-corrected exchange functional", Mol. Phys.
C>     <b>89</b> (1996) 433-445, DOI:
C>     <a href="https://doi.org/10.1080/002689796173813">
C>     10.1080/002689796173813</a>.
C>
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#if defined(NWAD_PRINT)
      Subroutine nwxc_x_gill_p(tol_rho, ipol, nq, wght, rho, rgamma, 
     &                         func)
#else
      Subroutine nwxc_x_gill(tol_rho, ipol, nq, wght, rho, rgamma, 
     &                       func)
#endif
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
      Subroutine nwxc_x_gill_d2(tol_rho, ipol, nq, wght, rho, rgamma, 
     &                          func)
#else
      Subroutine nwxc_x_gill_d3(tol_rho, ipol, nq, wght, rho, rgamma, 
     &                          func)
#endif
c
C$Id$
c
#include "nwad.fh"
c
      implicit none
c      
#include "nwxc_param.fh"
c
c     Input and other parameters
c
      double precision tol_rho !< [Input] The lower limit on the density
      integer ipol             !< [Input] The number of spin channels
      integer nq               !< [Input] The number of points
      double precision wght    !< [Input] The weight of the functional
c
c     Charge Density
c
      type(nwad_dble)::rho(nq,*)    !< [Input] The density
c
c     Charge Density Gradient
c
      type(nwad_dble)::rgamma(nq,*) !< [Input] The norm of the density gradients
c
c     Sampling Matrices for the XC Potential
c
      type(nwad_dble)::func(nq)    !< [Output] The value of the functional
c     double precision Amat(nq,*)   !< [Output] The derivative wrt rho
c     double precision Cmat(nq,*)   !< [Output] The derivative wrt rgamma
#ifdef SECOND_DERIV
c
c     Sampling Matrices for the XC Kernel
c
c     double precision Amat2(nq,*)  !< [Output] The 2nd derivative wrt rho
c     double precision Cmat2(nq,*)  !< [Output] The 2nd derivative wrt rgamma
c                                   !< and possibly rho
#endif
c
      double precision BETA, ONE3, FOUR3
#ifdef SECOND_DERIV
c     double precision SEVEN3
#endif
      Parameter (BETA = 1d0/137d0)
      Parameter (ONE3 = 1d0/3d0, FOUR3 = 4d0/3d0)
#ifdef SECOND_DERIV
c     Parameter (SEVEN3 = 7d0/3d0)
#endif
c
c References:
c
c    Gill , Mol. Phys. 89, 433 (1996)
c
c***************************************************************************
c
      integer n
      double precision C, d1x(2), d1g
      type(nwad_dble)::rhoval,rho13,rho43,gamma,x,g
#ifdef SECOND_DERIV
c     double precision d2x(3), d2g
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
            if (rho(n,R_T).lt.2.0d0*tol_rho) goto 10
            rhoval = 0.5d0*rho(n,R_T)
c
c           Spin alpha:
c
            rho13 = rhoval**ONE3
            rho43 = rho13*rhoval
c     Include factor of 4/3 in rho13 since it always appears with it
            rho13 = FOUR3*rho13
c
c           if (lfac) then
c              Ex = Ex + 2d0*rho43*C*qwght(n)*wght
               func(n) = func(n) + 2.d0*rho43*C*wght
c              Amat(n,D1_RA) = Amat(n,D1_RA) + rho13*C*wght
#ifdef SECOND_DERIV
c              Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA)
c    &              + ONE3*(rho13/rhoval)*C*wght
#endif
c           endif
c
            gamma = rgamma(n,G_TT)
c           gamma = delrho(n,1,1)*delrho(n,1,1) +
c    &              delrho(n,2,1)*delrho(n,2,1) +
c    &              delrho(n,3,1)*delrho(n,3,1)
            if (gamma.gt.tol_rho*tol_rho)then
               gamma = 0.25d0*gamma
               x = sqrt(gamma)/rho43
c              d1x(1) = -FOUR3*x/rhoval
c              d1x(2) = 0.5d0*x/gamma
               g = -BETA*x*sqrt(x)
c              d1g = -1.5d0*BETA*sqrt(x)
c
c              Ex = Ex + 2d0*rho43*g*qwght(n)*wght
               func(n) = func(n) + 2.d0*rho43*g*wght
c              Amat(n,D1_RA) = Amat(n,D1_RA)
c    &                       + (rho13*g+rho43*d1g*d1x(1))*wght
c              Cmat(n,D1_GAA) = Cmat(n,D1_GAA) + rho43*d1g*d1x(2)*wght
#ifdef SECOND_DERIV
c              d2g = 0.5d0*d1g/x
c              d2x(1) = -SEVEN3*d1x(1)/rhoval
c              d2x(2) = -FOUR3*d1x(2)/rhoval
c              d2x(3) = -0.5d0*d1x(2)/gamma
c              Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA)
c    &              + ONE3*(rho13/rhoval)*g*wght
c    &              + 2.d0*rho13*d1g*d1x(1)*wght
c    &              + rho43*(d2g*d1x(1)*d1x(1)+d1g*d2x(1))*wght
c              Cmat2(n,D2_RA_GAA) = Cmat2(n,D2_RA_GAA)
c    &              + rho13*d1g*d1x(2)*wght
c    &              + rho43*(d2g*d1x(1)*d1x(2)+d1g*d2x(2))*wght
c              Cmat2(n,D2_GAA_GAA) = Cmat2(n,D2_GAA_GAA)
c    &              + rho43*(d2g*d1x(2)*d1x(2)+d1g*d2x(3))*wght
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
            if (rho(n,R_A)+rho(n,R_B).lt.tol_rho) goto 20
            if (rho(n,R_A).lt.tol_rho) goto 25
c
c           Spin alpha:
c
            rhoval = rho(n,R_A)
            rho13 = rhoval**ONE3
            rho43 = rho13*rhoval
c     Include factor of 4/3 in rho13 since it always appears with it
            rho13 = FOUR3*rho13
c
c           if (lfac) then
c              Ex = Ex + rho43*C*qwght(n)*wght
               func(n) = func(n) + rho43*C*wght
c              Amat(n,D1_RA) = Amat(n,D1_RA) + rho13*C*wght
#ifdef SECOND_DERIV
c              Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA)
c    &              + ONE3*(rho13/rhoval)*C*wght
#endif
c           endif
c
            gamma = rgamma(n,G_AA)
c           gamma = delrho(n,1,1)*delrho(n,1,1) +
c    &              delrho(n,2,1)*delrho(n,2,1) +
c    &              delrho(n,3,1)*delrho(n,3,1)
            if (gamma.gt.tol_rho*tol_rho)then
               x = sqrt(gamma)/rho43
c              d1x(1) = -FOUR3*x/rhoval
c              d1x(2) = 0.5d0*x/gamma
               g = -BETA*x*sqrt(x)
c              d1g = -1.5d0*BETA*sqrt(x)
c
c              Ex = Ex + rho43*g*qwght(n)*wght
               func(n) = func(n) + rho43*g*wght
c              Amat(n,D1_RA) = Amat(n,D1_RA)
c    &                       + (rho13*g+rho43*d1g*d1x(1))*wght
c              Cmat(n,D1_GAA) = Cmat(n,D1_GAA) + rho43*d1g*d1x(2)*wght
#ifdef SECOND_DERIV
c              d2g = 0.5d0*d1g/x
c              d2x(1) = -SEVEN3*d1x(1)/rhoval
c              d2x(2) = -FOUR3*d1x(2)/rhoval
c              d2x(3) = -0.5d0*d1x(2)/gamma
c              Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA)
c    &              + ONE3*(rho13/rhoval)*g*wght
c    &              + 2.d0*rho13*d1g*d1x(1)*wght
c    &              + rho43*(d2g*d1x(1)*d1x(1)+d1g*d2x(1))*wght
c              Cmat2(n,D2_RA_GAA) = Cmat2(n,D2_RA_GAA)
c    &              + rho13*d1g*d1x(2)*wght
c    &              + rho43*(d2g*d1x(1)*d1x(2)+d1g*d2x(2))*wght
c              Cmat2(n,D2_GAA_GAA) = Cmat2(n,D2_GAA_GAA)
c    &              + rho43*(d2g*d1x(2)*d1x(2)+d1g*d2x(3))*wght
#endif
            endif
c
c
 25         continue
c
c           Spin beta:
c
            if (rho(n,R_B).lt.tol_rho) goto 20
c
            rhoval = rho(n,R_B)
            rho13 = rhoval**ONE3
            rho43 = rho13*rhoval
c     Include factor of 4/3 in rho13 since it always appears with it
            rho13 = FOUR3*rho13
c
c           if (lfac) then
c              Ex = Ex + rho43*C*qwght(n)*wght
               func(n) = func(n) + rho43*C*wght
c              Amat(n,D1_RB) = Amat(n,D1_RB) + rho13*C*wght
#ifdef SECOND_DERIV
c              Amat2(n,D2_RB_RB) = Amat2(n,D2_RB_RB)
c    &              + ONE3*(rho13/rhoval)*C*wght
#endif
c           endif
c
            gamma = rgamma(n,G_BB)
c           gamma = delrho(n,1,2)*delrho(n,1,2) +
c    &              delrho(n,2,2)*delrho(n,2,2) +
c    &              delrho(n,3,2)*delrho(n,3,2)
            if (gamma.gt.tol_rho*tol_rho)then
               x = sqrt(gamma)/rho43
c              d1x(1) = -FOUR3*x/rhoval
c              d1x(2) = 0.5d0*x/gamma
               g = -BETA*x*sqrt(x)
c              d1g = -1.5d0*BETA*sqrt(x)
c
c              Ex = Ex + rho43*g*qwght(n)*wght
               func(n) = func(n) + rho43*g*wght
c              Amat(n,D1_RB) = Amat(n,D1_RB)
c    &                       + (rho13*g+rho43*d1g*d1x(1))*wght
c              Cmat(n,D1_GBB) = Cmat(n,D1_GBB) + rho43*d1g*d1x(2)*wght
#ifdef SECOND_DERIV
c              d2g = 0.5d0*d1g/x
c              d2x(1) = -SEVEN3*d1x(1)/rhoval
c              d2x(2) = -FOUR3*d1x(2)/rhoval
c              d2x(3) = -0.5d0*d1x(2)/gamma
c              Amat2(n,D2_RB_RB) = Amat2(n,D2_RB_RB)
c    &              + ONE3*(rho13/rhoval)*g*wght
c    &              + 2.d0*rho13*d1g*d1x(1)*wght
c    &              + rho43*(d2g*d1x(1)*d1x(1)+d1g*d2x(1))*wght
c              Cmat2(n,D2_RB_GBB) = Cmat2(n,D2_RB_GBB)
c    &              + rho13*d1g*d1x(2)*wght
c    &              + rho43*(d2g*d1x(1)*d1x(2)+d1g*d2x(2))*wght
c              Cmat2(n,D2_GBB_GBB) = Cmat2(n,D2_GBB_GBB)
c    &              + rho43*(d2g*d1x(2)*d1x(2)+d1g*d2x(3))*wght
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
#ifndef NWAD_PRINT
#define NWAD_PRINT
c
c     Compile source again for Maxima
c
#include "nwxc_x_gill.F"
#endif
#ifndef SECOND_DERIV
#define SECOND_DERIV
c
c     Compile source again for the 2nd derivative case
c
#include "nwxc_x_gill.F"
#endif
#ifndef THIRD_DERIV
#define THIRD_DERIV
c
c     Compile source again for the 3rd derivative case
c
#include "nwxc_x_gill.F"
#endif
#undef NWAD_PRINT
C>
C> @}
