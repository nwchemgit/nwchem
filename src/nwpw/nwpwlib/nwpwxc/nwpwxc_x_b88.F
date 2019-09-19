#ifndef SECOND_DERIV
C> \ingroup nwpwxc
C> @{
C>
C> \file nwpwxc_x_b88.F
C> The B88 exchange functional
C>
C> @}
#endif
C> \ingroup nwpwxc_priv
C> @{
C>
C> \brief Evaluate the B88 exchange functional
C>
C> Evaluate the non-local part of B88 hybrid-GGA [1]. 
C>
C> ### References ###
C>
C> [1] A.D. Becke,
C> "Density-functional exchange-energy approximation with correct
C> asymptotic behavior",
C> Phys. Rev. A <b>38</b>, 3098-3100 (1998), DOI:
C> <a href="https://doi.org/10.1103/PhysRevA.38.3098">
C> 10.1103/PhysRevA.38.3098</a>.
C>
#if defined(FUJITSU_VPP)
!ocl scalar
#endif
#ifndef SECOND_DERIV
      Subroutine nwpwxc_x_b88(tol_rho, ipol, nq, wght, rho, rgamma,
     &                      func, Amat, Cmat)
#else
      Subroutine nwpwxc_x_b88_d2(tol_rho, ipol, nq, wght, rho,
     &                         rgamma, func, Amat, Amat2, Cmat, Cmat2)
#endif
c
C$Id$
c
      implicit none
c      
#include "nwpwxc_param.fh"
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
      double precision rho(nq,*)    !< [Input] The density
c
c     Charge Density Gradient
c
      double precision rgamma(nq,*) !< [Input] The norm of the density gradients
c
c     Sampling Matrices for the XC Potential
c
      double precision func(nq)    !< [Output] The value of the functional
      double precision Amat(nq,*)   !< [Output] The derivative wrt rho
      double precision Cmat(nq,*)   !< [Output] The derivative wrt rgamma
#ifdef SECOND_DERIV
c
c     Sampling Matrices for the XC Kernel
c
      double precision Amat2(nq,*)  !< [Output] The 2nd derivative wrt rho
      double precision Cmat2(nq,*)  !< [Output] The 2nd derivative wrt rgamma
                                    !< and possibly rho
#endif
c
      double precision BETA
      Parameter (BETA = 0.0042D0)

      double precision tolrho
      Parameter (tolrho=2.0d-8)
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
     &     dgdenom, t, ex,fx
#ifdef SECOND_DERIV
      double precision rhom23, d2g, d2gdenom
#endif
c
      arcsinh(x)=log(x+dsqrt(1d0+x*x))
      darcsinh(x)=1d0/dsqrt(1d0+x*x)
c
c     Uniform electron gas constant
c
      !lda_c =(3.0d0/2.0d0)*(3.0d0/(4.0d0*pi))**(1.0d0/3.0d0)
      C = -(1.5d0)*(0.75d0/acos(-1d0))**(1d0/3d0)
c
      if (ipol.eq.1) then
c
c        ======> SPIN-RESTRICTED <======
c
         do 10 n = 1, nq
            if (rho(n,R_T).lt.tol_rho) goto 10
c
c           Spin alpha:
c
            rho13 = (0.5d0*rho(n,R_T))**(1.d0/3.d0)
            rho43 = rho13**4 
c           gamma = delrho(n,1,1)*delrho(n,1,1) +
c    &              delrho(n,2,1)*delrho(n,2,1) +
c    &              delrho(n,3,1)*delrho(n,3,1)
            gamma = rgamma(n,G_TT)
            if (dsqrt(gamma).gt.tol_rho)then
               gamma = 0.25d0 * gamma
               x = dsqrt(gamma) / rho43
            else
               x = 0.0d0
            endif
c
            gdenom = 1d0 + 6d0*BETA*x*arcsinh(x)
            g = -BETA*x*x / gdenom
            dgdenom = 6d0*BETA*(arcsinh(x) + x*darcsinh(x))
            dg = BETA*x*(x*dgdenom - 2d0*gdenom) / gdenom**2

            ex = rho13*(-C + g)
            fx = (4d0/3d0)*rho13*(-C + (g-x*dg))
            if ((fx-ex).gt.0.0d0) then
               func(n) = func(n) + 2.d0*rho43*g*wght
               Amat(n,D1_RA) = Amat(n,D1_RA)
     &                       + (4d0/3d0)*rho13*(g-x*dg)*wght
           
               t = 0.5d0 * dg / sqrt(gamma) * wght
               Cmat(n,D1_GAA) = Cmat(n,D1_GAA) + t
c            Cmat(n,1,1) = Cmat(n,1,1) + t*delrho(n,1,1)
c            Cmat(n,2,1) = Cmat(n,2,1) + t*delrho(n,2,1)
c            Cmat(n,3,1) = Cmat(n,3,1) + t*delrho(n,3,1)
               
c
#ifdef SECOND_DERIV
c           if(lfac) g = g + C           ! Add local contribution back to g
               rhom23 = rho13 / (0.5d0*rho(n,R_T))
               d2gdenom = 6d0*BETA*darcsinh(x)*(2d0 - x*x/(x*x+1d0))
               d2g = -2d0*BETA/gdenom + 4d0*BETA*x*dgdenom/gdenom**2
     &              + BETA*x*x*d2gdenom/gdenom**2
     &              - 2d0*BETA*x*x*(dgdenom)**2/gdenom**3
c
               Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA)
     &              + (4d0/9d0)*rhom23*(g-x*dg+4d0*x*x*d2g)*wght
               Cmat2(n,D2_RA_GAA) = Cmat2(n,D2_RA_GAA)
     &              - (4d0/3d0)*(rhom23**2/rho(n,R_T))*d2g*wght
               Cmat2(n,D2_GAA_GAA) = Cmat2(n,D2_GAA_GAA)
     &              - 0.25d0*gamma**(-1.5d0)*(dg-x*d2g)*wght
#endif
            end if
c
 10      continue
c
      else
c
c        ======> SPIN-UNRESTRICTED <======
c
         do 20 n = 1, nq
c           if (rho(n,1).lt.tol_rho) goto 20
            if (rho(n,R_A).lt.0.5d0*tol_rho) goto 25
c
c           Spin alpha:
c
            rho13 = rho(n,R_A)**(1.d0/3.d0)
            rho43 = rho13*rho(n,R_A)
c           gamma = delrho(n,1,1)*delrho(n,1,1) +
c    &              delrho(n,2,1)*delrho(n,2,1) +
c    &              delrho(n,3,1)*delrho(n,3,1)
            gamma = rgamma(n,G_AA)
            if (dsqrt(gamma).gt.tol_rho)then
               x = dsqrt(gamma) / rho43
            else
               x = 0d0
            endif
c
            gdenom = 1d0 + 6d0*BETA*x*arcsinh(x)
            g = -BETA*x*x / gdenom
            dgdenom = 6d0*BETA*(arcsinh(x) + x*darcsinh(x))
            dg = BETA*x*(x*dgdenom - 2d0*gdenom) / gdenom**2
c
            ex = rho13*(-C + g)
            fx = (4d0/3d0)*rho13*(-C + (g-x*dg))
            if ((fx-ex).gt.0.0d0) then
               func(n) = func(n) + rho43*g*wght
               Amat(n,D1_RA) = Amat(n,D1_RA) 
     &                       + (4d0/3d0)*rho13*(g-x*dg)*wght
c
               t = dg / sqrt(gamma) * wght
               Cmat(n,D1_GAA) = Cmat(n,D1_GAA) + t * 0.5d0
c               Cmat(n,1,1) = Cmat(n,1,1) + t*delrho(n,1,1)
c               Cmat(n,2,1) = Cmat(n,2,1) + t*delrho(n,2,1)
c               Cmat(n,3,1) = Cmat(n,3,1) + t*delrho(n,3,1)
c
#ifdef SECOND_DERIV
               rhom23 = rho13 / rho(n,R_A)
               d2gdenom = 6d0*BETA*darcsinh(x)*(2d0 - x*x/(x*x+1d0))
               d2g = -2d0*BETA/gdenom + 4d0*BETA*x*dgdenom/gdenom**2
     &           + BETA*x*x*d2gdenom/gdenom**2
     &           - 2d0*BETA*x*x*(dgdenom)**2/gdenom**3
c
               Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA)
     &           + (4d0/9d0)*rhom23*(g-x*dg+4d0*x*x*d2g)*wght
               Cmat2(n,D2_RA_GAA) = Cmat2(n,D2_RA_GAA)
     &           - (2d0/3d0)*(rhom23**2/rho(n,R_A))*d2g*wght
               if (x.gt.tol_rho) then
               Cmat2(n,D2_GAA_GAA) = Cmat2(n,D2_GAA_GAA)
     &              - 0.25d0*gamma**(-1.5d0)*(dg-x*d2g)*wght
               endif
#endif
            endif
c
 25         continue
c
c           Spin beta:
c
            if (rho(n,R_B).lt.0.5d0*tol_rho) goto 20
c
            rho13 = rho(n,R_B)**(1.d0/3.d0)
            rho43 = rho13*rho(n,R_B)
c           gamma = delrho(n,1,2)*delrho(n,1,2) +
c    &              delrho(n,2,2)*delrho(n,2,2) +
c    &              delrho(n,3,2)*delrho(n,3,2)
            gamma = rgamma(n,G_BB)
            if (dsqrt(gamma).gt.tol_rho)then
               x = dsqrt(gamma) / rho43
            else
               x = 0d0
            endif
c
            gdenom = 1d0 + 6d0*BETA*x*arcsinh(x)
            g = -BETA*x*x / gdenom
            dgdenom = 6d0*BETA*(arcsinh(x) + x*darcsinh(x))
            dg = BETA*x*(x*dgdenom - 2d0*gdenom) / gdenom**2
c
            ex = rho13*(-C + g)
            fx = (4d0/3d0)*rho13*(-C + (g-x*dg))
            if ((fx-ex).gt.0.0d0) then
               func(n) = func(n) +rho43*g*wght
               Amat(n,D1_RB) = Amat(n,D1_RB)
     &                       + (4d0/3d0)*rho13*(g-x*dg)*wght
c
               t = dg / sqrt(gamma) * wght
               Cmat(n,D1_GBB) = Cmat(n,D1_GBB) + t * 0.5d0
c               Cmat(n,1,2) = Cmat(n,1,2) + t*delrho(n,1,2)
c               Cmat(n,2,2) = Cmat(n,2,2) + t*delrho(n,2,2)
c               Cmat(n,3,2) = Cmat(n,3,2) + t*delrho(n,3,2)
c
#ifdef SECOND_DERIV
               rhom23 = rho13 / rho(n,R_B)
               d2gdenom = 6d0*BETA*darcsinh(x)*(2d0 - x*x/(x*x+1d0))
               d2g = -2d0*BETA/gdenom + 4d0*BETA*x*dgdenom/gdenom**2
     &           + BETA*x*x*d2gdenom/gdenom**2
     &           - 2d0*BETA*x*x*(dgdenom)**2/gdenom**3
c
               Amat2(n,D2_RB_RB) = Amat2(n,D2_RB_RB)
     &           + (4d0/9d0)*rhom23*(g-x*dg+4d0*x*x*d2g)*wght
               Cmat2(n,D2_RB_GBB) = Cmat2(n,D2_RB_GBB)
     &           - (2d0/3d0)*(rhom23**2/rho(n,R_B))*d2g*wght
               if (x.gt.tol_rho) then
               Cmat2(n,D2_GBB_GBB) = Cmat2(n,D2_GBB_GBB)
     &              - 0.25d0*gamma**(-1.5d0)*(dg-x*d2g)*wght
               endif
#endif
            endif
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
#include "nwpwxc_x_b88.F"
#endif
C>
C> @}
