#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#if !defined(NWAD_PRINT)
C> \ingroup nwxc
C> @{
C>
C> \file nwxc_xc_kt1.F
C> The KT1, KT2 and SSB-D functional
C>
C> @}
#endif
#endif
C>
C> \ingroup nwxc_priv
C> @{
C>
C> \brief Evaluate the KT1, KT2, and part of the SSB-D functional
C>
C> This is the gradient correction term by Keal and Tozer [1],
C> which can be used on its own (as in the KT1 or KT2 functional),
C> or forms part of the SSB-D functional [2].
C>
C> Note that even though the energy takes the form of an exchange
C> functional, this gradient correction is strictly speaking
C> NOT an exchange functional!
C>
C> ### References ###
C>
C> [1] T.W. Keal, D.J. Tozer,
C>     "The exchange-correlation potential in Kohn-Sham nuclear
C>     magnetic resonance shielding calculations", J. Chem. Phys.
C>     <b>119</b>, 3015-3024 (2006), DOI:
C>     <a href="https://doi.org/10.1063/1.1590634">
C>     10.1063/1.1590634</a>.
C>
C> [2] M. Swart, M. Sola, F.M. Bickelhaupt,
C>     "A new all-round density functional based on spin states and
C>     SN2 barriers",
C>     J. Chem. Phys. <b>131</b>, 094103 (2009), DOI:
C>     <a href="https://doi.org/10.1063/1.3213193">
C>     10.1063/1.3213193</a>.
C>
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#if defined(NWAD_PRINT)
      Subroutine nwxc_xc_kt1_p(tol_rho, ipol, nq, wght, rho, rgamma, 
     &                         func)
#else
      Subroutine nwxc_xc_kt1(tol_rho, ipol, nq, wght, rho, rgamma, 
     &                       func)
#endif
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
      Subroutine nwxc_xc_kt1_d2(tol_rho, ipol, nq, wght, rho, rgamma, 
     &                          func)
#else
      Subroutine nwxc_xc_kt1_d3(tol_rho, ipol, nq, wght, rho, rgamma, 
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
      type(nwad_dble)::func(nq)     !< [Output] The value of the functional
c     double precision Amat(nq,*)   !< [Output] The derivative wrt rho
c     double precision Cmat(nq,*)   !< [Output] The derivative wrt rgamma

#ifdef SECOND_DERIV
c
c     Second Derivatives of the Exchange Energy Functional
c
c     double precision Amat2(nq,*)  !< [Output] The 2nd derivative wrt rho
c     double precision Cmat2(nq,*)  !< [Output] The 2nd derivative wrt rgamma
c                                   !< and possibly rho
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
      type(nwad_dble)::hrho, rho13, rho43, gamma
      type(nwad_dble)::g, gdenom, gdenom2
#ifdef SECOND_DERIV
c     double precision rho23, rhom23, gdenom3
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
c           if (sqrt(gamma).gt.tol_rho) then
               gamma = 0.25d0 * gamma
c           else
c              goto 10
c           endif
c
            gdenom = 1.d0 / (rho43 + DELTA)
            gdenom2 = gdenom*gdenom
            g = GAMKT * gamma * gdenom
c
            func(n) = func(n) + 2.d0*g*wght
c           Amat(n,D1_RA) = Amat(n,D1_RA) - (4d0/3d0)*GAMKT*gamma*rho13*
c    &                  wght*gdenom2
c           Cmat(n,D1_GAA) = Cmat(n,D1_GAA) + GAMKT*gdenom*wght
c
#ifdef SECOND_DERIV
c           rho23 = rho13*rho13
c           rhom23 = rho13 / (0.5d0*rho(n,R_T))
c           gdenom3 = gdenom2*gdenom
c
c           Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA) + (4d0/3d0)*GAMKT*
c    &          gamma*(rho23*7d0/3d0 - DELTA*rhom23/3d0)*gdenom3*wght
c           Cmat2(n,D2_RA_GAA) = Cmat2(n,D2_RA_GAA)
c    &           - (4d0/3d0)*GAMKT*rho13*gdenom2*wght
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
c           if (sqrt(gamma).lt.tol_rho) then
c              goto 25
c           endif
c
            gdenom = 1d0 / (rho43 + DELTA)
            gdenom2 = gdenom*gdenom
            g = GAMKT * gamma * gdenom
c
            func(n) = func(n) + g*wght
c           Amat(n,D1_RA) = Amat(n,D1_RA) - (4d0/3d0)*GAMKT*gamma*rho13*
c    &                 gdenom2*wght
c           Cmat(n,D1_GAA) = Cmat(n,D1_GAA) + GAMKT*gdenom*wght
c
#ifdef SECOND_DERIV
c           rho23 = rho13*rho13
c           rhom23 = rho13 / rho(n,R_A)
c           gdenom3 = gdenom2*gdenom
c
c           Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA) + (4d0/3d0)*GAMKT*
c    &          gamma*(rho23*7d0/3d0 - DELTA*rhom23/3d0)*gdenom3*wght
c           Cmat2(n,D2_RA_GAA) = Cmat2(n,D2_RA_GAA)
c    &           - (4d0/3d0)*GAMKT*rho13*gdenom2*wght
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
            if (rho(n,R_B).lt.tol_rho) goto 20
c
            rho13 = rho(n,R_B)**(1.d0/3.d0)
            rho43 = rho13*rho(n,R_B)
            gamma = rgamma(n,G_BB)
c           gamma = delrho(n,1,2)*delrho(n,1,2) +
c    &              delrho(n,2,2)*delrho(n,2,2) +
c    &              delrho(n,3,2)*delrho(n,3,2)
c           if (sqrt(gamma).lt.tol_rho) then
c              goto 20
c           endif
c           
            gdenom = 1d0 / (rho43 + DELTA)
            gdenom2 = gdenom*gdenom 
            g = GAMKT * gamma * gdenom
c           
            func(n) = func(n) + g*wght
c           Amat(n,D1_RB) = Amat(n,D1_RB) - (4d0/3d0)*GAMKT*gamma*rho13*
c    &            wght*gdenom2
c           Cmat(n,D1_GBB) = Cmat(n,D1_GBB) + GAMKT*gdenom*wght
c           
#ifdef SECOND_DERIV
c           rho23 = rho13*rho13
c           rhom23 = rho13 / rho(n,R_B)
c           gdenom3 = gdenom2*gdenom
c
c           Amat2(n,D2_RB_RB) = Amat2(n,D2_RB_RB) + (4d0/3d0)*GAMKT*
c    &          gamma*(rho23*7d0/3d0 - DELTA*rhom23/3d0)*gdenom3*wght
c           Cmat2(n,D2_RB_GBB) = Cmat2(n,D2_RB_GBB)
c    &           - (4d0/3d0)*GAMKT*rho13*gdenom2*wght
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
#ifndef NWAD_PRINT
#define NWAD_PRINT
c
c     Compile source again for the 2nd derivative case
c
#include "nwxc_xc_kt1.F"
#endif
#ifndef SECOND_DERIV
#define SECOND_DERIV
c
c     Compile source again for the 2nd derivative case
c
#include "nwxc_xc_kt1.F"
#endif
#ifndef THIRD_DERIV
#define THIRD_DERIV
c
c     Compile source again for the 3rd derivative case
c
#include "nwxc_xc_kt1.F"
#endif
#undef NWAD_PRINT
C>
C> @}
