C> \ingroup nwxc
C> @{
C>
C> \file nwxc_cs00.F
C> The Casida-Salahub potential correction
C>
C> @}
C>
C> \ingroup nwxc_priv
C> @{
C>
C> \brief Evaluate the Casida-Salahub potential correction
C>
C> Evaluate the Casida-Salahub potential correction [1]. This correction
C> specifically targets the long range interaction and improves
C> Rydberg states.
C>
C> ### References ###
C>
C> [1] M.E. Casida, D.R. Salahub,
C>     "Asymptotic correction approach to improving approximate
C>     exchange–correlation potentials: Time-dependent
C>     density-functional theory calculations of molecular excitation
C>     spectra", J. Chem. Phys. <b>113</b>, 8918-8935 (2000), DOI:
C>     <A HREF="https://doi.org/10.1063/1.1319649">
C>     10.1063/1.1319649</A>.
C>
#if defined(FUJITSU_VPP)
!ocl scalar
#endif
      Subroutine nwxc_cs00(param, tol_rho, ipol, nq, wght, rho,
     &                     rgamma, func, Amat, Cmat)
c
C$Id$
c
      implicit none
c      
#include "nwxc_param.fh"
c
      double precision param(*)!< [Input] Parameters of the functional
                               !< - param(1): delta
                               !< - param(2): \f$ \epsilon_{\mathrm{homo}} \f$
      double precision tol_rho !< [Input] The lower limit on the density
      integer nq               !< [Input] The number of points
      integer ipol             !< [Input] The number of spin channels
      double precision wght    !< [Input] The weight of the functional
c
c     Charge Density
c
      double precision rho(nq,*) !< [Input] The density
c
c     Charge Density Gradient
c
      double precision rgamma(nq,*) !< [Input] The norm of the density gradients
c
c     Sampling Matrices for the XC Potential & Energy
c
      double precision func(nq)   !< [Output] The value of the functional
      double precision Amat(nq,*) !< [Output] The derivative wrt rho
      double precision Cmat(nq,*) !< [Output] The derivative wrt rgamma
c
      double precision delta, shift, e_homo
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
      arcsinh(x)=log(x+dsqrt(1d0+x*x))
c
      delta  = param(1)
      e_homo = param(2)
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
            if (rho(n,R_T).lt.tol_rho) goto 10
c
c           Spin alpha:
c
            hrho  = 0.5d0*rho(n,R_T)
            rho13 = hrho**(1.d0/3.d0)
            rho43 = rho13*hrho
            gamma = rgamma(n,G_TT)
c           gamma = delrho(n,1,1)*delrho(n,1,1) +
c    &              delrho(n,2,1)*delrho(n,2,1) +
c    &              delrho(n,3,1)*delrho(n,3,1)
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
            p1 = Amat(n,R_A) + rho13*g*(1d0-wght)
            p2 = Amat(n,R_A) - shift
            Amat(n,R_A) = max(p1, p2)
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
            rho13 = rho(n,R_A)**(1.d0/3.d0)
            rho43 = rho13*rho(n,R_A)
            gamma = rgamma(n,G_AA)
c           gamma = delrho(n,1,1)*delrho(n,1,1) +
c    &              delrho(n,2,1)*delrho(n,2,1) +
c    &              delrho(n,3,1)*delrho(n,3,1)
            if (dsqrt(gamma).gt.tol_rho)then
               x = dsqrt(gamma) / rho43
            else
               x = 0d0
            endif
c
            gdenom = 1d0 + 3d0*BETA*x*arcsinh(x)
            g = -BETA*x*x / gdenom
c
            p1 = Amat(n,R_A) + rho13*g*(1d0-wght)
            p2 = Amat(n,R_A) - shift
            Amat(n,R_A) = max(p1, p2)
c
 25         continue
c
c           Spin beta:
c
            if (rho(n,R_B).lt.tol_rho) goto 20
c
            rho13 = rho(n,R_B)**(1.d0/3.d0)
            rho43 = rho13*rho(n,R_B)
            gamma = rgamma(n,G_BB)
c           gamma = delrho(n,1,2)*delrho(n,1,2) +
c    &              delrho(n,2,2)*delrho(n,2,2) +
c    &              delrho(n,3,2)*delrho(n,3,2)
            if (dsqrt(gamma).gt.tol_rho)then
               x = dsqrt(gamma) / rho43
            else
               x = 0d0
            endif
c
            gdenom = 1d0 + 3d0*BETA*x*arcsinh(x)
            g = -BETA*x*x / gdenom
c
            p1 = Amat(n,R_B) + rho13*g*(1d0-wght)
            p2 = Amat(n,R_B) - shift
            Amat(n,R_B) = max(p1, p2)
c
 20      continue
c
      endif
c
      return
      end
C> @}
