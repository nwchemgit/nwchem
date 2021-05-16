#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#if !defined(NWAD_PRINT)
C> \ingroup nwxc
C> @{
C>
C> \file nwxc_c_perdew86.F
C> The Perdew correlation functional of 1986
C>
C> @}
#endif
#endif
C> \ingroup nwxc_priv
C> @{
C>
C> \brief Evaluate the Perdew 1986 correlation functional
C>
C> Evaluates the Perdew 1986 GGA correlation functional [1,2,3].
C>
C> ### References ###
C>
C> [1] J.P. Perdew,
C>     "Density-functional approximation for the correlation energy of
C>     the inhomogeneous electron gas", Phys. Rev. B <b>33</b>,
C>     8822–8824 (1986), DOI:
C>     <a href="https://doi.org/10.1103/PhysRevB.33.8822">
C>     10.1103/PhysRevB.33.8822</a>.
C>
C> [2] P. Mlynarski, D.R. Salahub,
C>     "Self-consistent implementation of nonlocal exchange and
C>     correlation in a Gaussian density-functional method",
C>     Phys. Rev. B <b>43</b>, 1399–1410 (1991), DOI:
C>     <a href="https://doi.org/10.1103/PhysRevB.43.1399">
C>     10.1103/PhysRevB.43.1399</a>.
C>
C> [3] J.P. Perdew,
C>     "Erratum: Density-functional approximation for the correlation
C>     energy of the inhomogeneous electron gas", Phys. Rev. B
C>     <b>34</b>, 7406–7406 (1986), DOI:
C>     <a href="https://doi.org/10.1103/PhysRevB.34.7406">
C>     10.1103/PhysRevB.34.7406</a>.
C>
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#if defined(NWAD_PRINT)
      Subroutine nwxc_c_perdew86_p(tol_rho, ipol, nq, wght, rho, rgamma,
     &                             ffunc)
#else
      Subroutine nwxc_c_perdew86(tol_rho, ipol, nq, wght, rho, rgamma, 
     &                           ffunc)
#endif
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
      Subroutine nwxc_c_perdew86_d2(tol_rho, ipol, nq, wght,
     &                              rho, rgamma, ffunc) 
#else
      Subroutine nwxc_c_perdew86_d3(tol_rho, ipol, nq, wght,
     &                              rho, rgamma, ffunc)
#endif
c
c$Id$
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
      type(nwad_dble)::ffunc(nq)    !< [Output] The value of the functional
c     double precision Amat(nq,*)   !< [Output] The derivative wrt rho
c     double precision Cmat(nq,*)   !< [Output] The derivative wrt rgamma
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
c
c     Sampling Matrices for the XC Kernel
c
c     double precision Amat2(nq,*)  !< [Output] The 2nd derivative wrt rho
c     double precision Cmat2(nq,*)  !< [Output] The 2nd derivative wrt rgamma
c                                   !< and possibly rho
#endif
#if defined(THIRD_DERIV)
c
c     Sampling Matrices for the XC Kernel
c
c     double precision Amat3(nq,*)  !< [Output] The 3rd derivative wrt rho
c     double precision Cmat3(nq,*)  !< [Output] The 3rd derivative wrt rgamma
c                                   !< and possibly rho
#endif
      double precision TOLL, EXPTOL, alpha, beta, pgamma, delta, 
     &                 beta10, ftilde, zzz, fff, pfff, CINF, ONE, 
     &                 ONE3, THREE, FOUR3, SEV6, FIVE3, 
     &                 TWO3, FIVE6, pi
      double precision SEVEN3, EIGHT3
      Parameter (TOLL = 1.D-40, EXPTOL = 80.d0)
      Parameter (alpha = 0.023266D0, beta  =  7.389D-6, 
     &   pgamma = 8.723d0, delta = 0.472d0,  beta10 = 10000.d0*beta)
      parameter (ftilde = 0.11d0, zzz = 0.001667d0, fff = 0.002568d0)
      parameter(pfff = 1.745d0, CINF = zzz+fff)
      Parameter (ONE = 1.D0, ONE3 = 1.d0/3.d0, THREE = 3.d0)
      Parameter (FOUR3 = 4.D0/3.D0, SEV6 = 7.d0/6.d0)
      parameter (FIVE3 = 5.d0/3.d0, TWO3 = 2.d0/3.d0, FIVE6 = 5.d0/6.d0)
      parameter (SEVEN3 = 7.0d0/3.0d0, EIGHT3 = 8.0d0/3.0d0)
c     parameter (pi = 3.1415926535897932385d0)
c     
c     Mlynarski Salahub PRB 43, 1399 (1991)
c     
      integer n
      type(nwad_dble)::rhoval
      double precision rsfact
      type(nwad_dble)::rs, rs2, rs3
      type(nwad_dble)::rho13, rho43, rho76, arho
      double precision d1rs
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
c     double precision d2rs
#endif
#if defined(THIRD_DERIV)
c     double precision d3rs
#endif
      type(nwad_dble)::gamma,gam12,zeta,func,dt12,phi,d,expfac
      type(nwad_dble)::anum,aden,Cn,dm1
      double precision d1anum, d1aden, d1Cn,
     &     d1phi(2), dlnphi, d1f(3),
     &     dlnfrho(2), dlnfgam
      double precision d1z(2), adp, d1d(2), t,
     &     d1dt12
      double precision aden2
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
c     double precision d2anum, d2aden, rrho2, d2z(3), dpp, d2d(3),
c    &     d2phi(3), d2dt12, d2Cn
c     double precision aden3
c     double precision arho2
c     double precision d2lnphi
c     double precision d2f(3)
c     double precision d2lnfrho(3), d2lnfrg(2), d2lnfgam
#endif
#if defined(THIRD_DERIV)
c     double precision d3lnphi
c     double precision d3anum, d3aden, d3Cn, d3phi(4)
c     double precision d3lnfrho(4), d3lnfgam
c     double precision d3f(3)
c     double precision aden4
c     double precision arho3
#endif
c
      pi = acos(-1.0d0)
      rsfact = (0.75d0/pi)**ONE3
c
      if (ipol.eq.1 )then
c
c        ======> SPIN-RESTRICTED <======
c
         do 10 n = 1, nq
            rhoval = rho(n,R_T)
            if (rhoval.lt.tol_rho) goto 10
            arho=1.d0/rhoval
            rho13 = rhoval**ONE3
            rho43 = rhoval*rho13
            rho76 = rhoval**SEV6
            rs = rsfact/rho13
            rs2 = rs*rs
            rs3 = rs2*rs
c           d1rs = -ONE3*rs*arho
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
c           d2rs = -FOUR3*d1rs*arho
#endif
#if defined(THIRD_DERIV)
c           d3rs = -SEVEN3*d2rs*arho 
#endif
            gamma = rgamma(n,G_TT)
c           gamma = delrho(n,1,1)*delrho(n,1,1) +
c    &              delrho(n,2,1)*delrho(n,2,1) +
c    &              delrho(n,3,1)*delrho(n,3,1)
            if (gamma.gt.tol_rho*tol_rho) then
              gam12 = sqrt(gamma)
            else
              gam12 = 0.0d0
            endif
c         
c           C(n)
c         
            anum = fff+alpha*rs+beta*rs2
            aden = 1.d0+pgamma*rs+delta*rs2+beta10*rs3
            Cn = zzz + anum/aden
c           d1anum = alpha + 2d0*beta*rs
c           d1aden = pgamma + 2d0*delta*rs + 3d0*beta10*rs2
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
c           d2anum = 2d0*beta
c           d2aden = 2d0*delta + 6d0*beta10*rs
#endif
#if defined(THIRD_DERIV)
c           d3anum = 0.0d0
c           d3aden = 6.0d0*beta10 
#endif
c     First compute rs derivative
c           aden2 = aden*aden
c           d1Cn = d1anum/aden - anum*d1aden/aden2
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
c           aden3 = aden2*aden
c           d2Cn = d2anum/aden - (2d0*d1anum*d1aden+anum*d2aden)/aden2
c    &           + 2d0*anum*d1aden**2/aden3
#endif
#if defined(THIRD_DERIV)
c           aden4 = aden3*aden
c
c           d3Cn = -( 3.0d0*d2anum*d1aden + 3.0d0*d1anum*d2aden
c    1              + anum*d3aden )/aden2
c    2           + 6.0d0*( d1anum*d1aden**2
c    3                   + anum*d2aden*d1aden )/aden3
c    4           - 6.0d0*anum*d1aden**3/aden4
#endif
c     Convert to rho derivative
#if defined(THIRD_DERIV)
c           d3Cn = d3Cn*d1rs*d1rs*d1rs
c    1           + 3.0d0*d2Cn*d2rs*d1rs
c    2           + d1Cn*d3rs
#endif
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
c           d2Cn = d2Cn*d1rs*d1rs + d1Cn*d2rs
#endif
c           d1Cn = d1Cn*d1rs
c         
c           phi(n,gradn)
c         
            expfac = 0.d0
            phi = (pfff*ftilde)*(CINF/Cn)*gam12/rho76
            if (phi.lt.EXPTOL) expfac = exp(-phi)
c           dlnphi = -(d1Cn/Cn + SEV6/rhoval)
c           d1phi(1) = phi*dlnphi
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
c           arho2 = arho*arho
c           d2lnphi = (d1Cn/Cn)**2 - d2Cn/Cn + SEV6*arho2
c           d2phi(1) = d1phi(1)*dlnphi + phi*d2lnphi
c           d2phi(1) = d1phi(1)*dlnphi
c    &               + phi*((d1Cn/Cn)**2 - d2Cn/Cn + SEV6/rhoval**2)
#endif
#if defined(THIRD_DERIV)
c           arho3 = arho2*arho
c
c           d3lnphi = -2.0d0*(d1Cn/Cn)**3
c    1              + 3.0d0*(d2Cn/Cn)*(d1Cn/Cn)
c    2              - d3Cn/Cn
c    3              - SEVEN3*arho3
c           d3phi(1) = d2phi(1)*dlnphi
c    1               + 2.0d0*d1phi(1)*d2lnphi
c    2               + phi*d3lnphi
#endif
c         
c           functional
c         
            func = expfac*Cn*gamma/rho43
c           dlnfrho(1) = d1Cn/Cn - (d1phi(1) + FOUR3/rhoval)
c           d1f(1) = dlnfrho(1)*func
c           Amat(n,D1_RA) = Amat(n,D1_RA) + d1f(1)*wght
c           if (gam12.gt.TOLL)then
c              d1phi(2) = phi / (2d0*gamma)
c              dlnfgam = 1d0/gamma - d1phi(2)
c              d1f(3) = func*dlnfgam
c              Cmat(n,D1_GAA) = Cmat(n,D1_GAA) + d1f(3)*wght
c              Cmat(n,D1_GAB) = Cmat(n,D1_GAB) + d1f(3)*2D0*wght
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
c              d2phi(2) = d1phi(2)*dlnphi
c              d2phi(3) =-d1phi(2)/(2d0*gamma)
c!!! Which of the following are actually needed for restricted?
c!!! Should treat derivatives of d as zero? d is a constant?
c Daniel (11-19-12): d is a constant (it equals 1) for a restricted 
c calculation, since there is no spin-polarization.  Thus, the
c derivatives are zero.
c              d2lnfrho(1) = -d2phi(1) - (d1Cn/Cn)**2 + d2Cn/Cn
c    1                     + FOUR3*arho2
c
c              d2f(1) = d1f(1)*dlnfrho(1)
c    1                + func*d2lnfrho(1)
c
c              t = d2f(1)*wght
c
c              Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA) + t
c              Amat2(n,D2_RA_RB) = Amat2(n,D2_RA_RB) + t
c    &              + (d1f(1)*dlnfrho(1)
c    &              + func*t)*wght
#if 0
c              Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA)
c    &              + (d1f(1)*dlnfrho(1)
c    &              + func*(d1d(1)*d1d(1)*dm1**2-d2d(1)*dm1+t))*wght
c              Amat2(n,D2_RA_RB) = Amat2(n,D2_RA_RB)
c    &              + (d1f(1)*dlnfrho(2)
c    &              + func*(d1d(1)*d1d(2)*dm1**2-d2d(2)*dm1+t))*wght
#endif
c rg terms
c              d2lnfrg(1) = -d2phi(2)
c              d2f(2) = (d1f(1)*dlnfgam + func*d2lnfrg(1))
c              t = d2f(2)*wght
c
c              Cmat2(n,D2_RA_GAA) = Cmat2(n,D2_RA_GAA) + t
c              Cmat2(n,D2_RA_GAB) = Cmat2(n,D2_RA_GAB) + t*2d0
c              Cmat2(n,D2_RA_GBB) = Cmat2(n,D2_RA_GBB) + t
c gg terms
c              d2lnfgam = -1.0d0/gamma**2 - d2phi(3)
c              d2f(3) = d1f(3)*dlnfgam + func*d2lnfgam
c              t = d2f(3)*wght
c
c              Cmat2(n,D2_GAA_GAA) = Cmat2(n,D2_GAA_GAA) + t
c              Cmat2(n,D2_GAA_GBB) = Cmat2(n,D2_GAA_GBB) + t
c              Cmat2(n,D2_GAA_GAB) = Cmat2(n,D2_GAA_GAB) + t*2d0
c              Cmat2(n,D2_GAB_GAB) = Cmat2(n,D2_GAB_GAB) + t*4d0
#endif
#if defined(THIRD_DERIV)
c rrr terms
c              d3lnfrho(1) = -d3phi(1)
c    1                     + 2.0d0*(d1Cn/Cn)**3
c    2                     - 3.0d0*(d2Cn/Cn)*(d1Cn/Cn)
c    3                     + d3Cn/Cn
c    4                     - EIGHT3*arho3
c
c              d3f(1) = d2f(1)*dlnfrho(1)
c    1                + 2.0d0*d1f(1)*d2lnfrho(1)
c    2                + func*d3lnfrho(1)
c
c              t = d3f(1)*wght
c
c              Amat3(n,D3_RA_RA_RA) = Amat3(n,D3_RA_RA_RA) + t
c              Amat3(n,D3_RA_RA_RB) = Amat3(n,D3_RA_RA_RB) + t
c              Amat3(n,D3_RA_RB_RB) = Amat3(n,D3_RA_RB_RB) + t
c rrg terms
c              d3phi(2) = d2phi(2)*dlnphi + d1phi(2)*d2lnphi
c
c              t = ( d2f(2)*dlnfrho(1)
c    1             - d1f(1)*d2phi(2)
c    2             + d1f(3)*d2lnfrho(1)
c    3             - func*d3phi(2) )*wght
c
c              Cmat3(n,D3_RA_RA_GAA) = Cmat3(n,D3_RA_RA_GAA) + t
c              Cmat3(n,D3_RA_RA_GAB) = Cmat3(n,D3_RA_RA_GAB) + t*2.0d0
c              Cmat3(n,D3_RA_RA_GBB) = Cmat3(n,D3_RA_RA_GBB) + t
c              Cmat3(n,D3_RA_RB_GAA) = Cmat3(n,D3_RA_RB_GAA) + t
c              Cmat3(n,D3_RA_RB_GAB) = Cmat3(n,D3_RA_RB_GAB) + t*2.0d0
c              Cmat3(n,D3_RA_RB_GBB) = Cmat3(n,D3_RA_RB_GBB) + t
c rgg terms
c              d3phi(3) = -d2phi(3)*dlnphi
c
c              t = ( d2f(2)*dlnfgam
c    1             + d1f(1)*d2lnfgam
c    2             + d1f(3)*d2lnfrg(1)
c    3             + func*d3phi(3) )*wght
c
c              Cmat3(n,D3_RA_GAA_GAA) = Cmat3(n,D3_RA_GAA_GAA) + t
c              Cmat3(n,D3_RA_GAA_GAB) = Cmat3(n,D3_RA_GAA_GAB) + t*2.0d0
c              Cmat3(n,D3_RA_GAA_GBB) = Cmat3(n,D3_RA_GAA_GBB) + t
c              Cmat3(n,D3_RA_GAB_GAB) = Cmat3(n,D3_RA_GAB_GAB) + t*4.0d0
c              Cmat3(n,D3_RA_GAB_GBB) = Cmat3(n,D3_RA_GAB_GBB) + t*2.0d0
c              Cmat3(n,D3_RA_GBB_GBB) = Cmat3(n,D3_RA_GBB_GBB) + t
c ggg terms              
c              d3phi(4) = -3.0d0*d2phi(3)/(2.0d0*gamma)
c              d3lnfgam = 2.0d0/gamma**3 - d3phi(4)
c
c              t = ( d2f(3)*dlnfgam
c    1             + 2.0d0*d1f(3)*d2lnfgam
c    2             + func*d3lnfgam )*wght
c
c              Cmat3(n,D3_GAA_GAA_GAA) = Cmat3(n,D3_GAA_GAA_GAA) + t
c              Cmat3(n,D3_GAA_GAA_GAB) = Cmat3(n,D3_GAA_GAA_GAB) 
c    1                                 + t*2.0d0
c              Cmat3(n,D3_GAA_GAA_GBB) = Cmat3(n,D3_GAA_GAA_GBB) + t
c              Cmat3(n,D3_GAA_GAB_GAB) = Cmat3(n,D3_GAA_GAB_GAB) 
c    1                                 + t*4.0d0
c              Cmat3(n,D3_GAA_GAB_GBB) = Cmat3(n,D3_GAA_GAB_GBB) 
c    1                                 + t*2.0d0
c              Cmat3(n,D3_GAA_GBB_GBB) = Cmat3(n,D3_GAA_GBB_GBB) + t
c              Cmat3(n,D3_GAB_GAB_GAB) = Cmat3(n,D3_GAB_GAB_GAB) 
c    1                                 + t*8.0d0
#endif
c           endif
            ffunc(n)=ffunc(n)+func*wght
   10    continue
      else
c
c        ======> SPIN-UNRESTRICTED <======
c
         do 20 n = 1, nq
            rhoval = 0.0d0
            gamma  = 0.0d0
            if (rho(n,R_A).ge.0.5d0*tol_rho) then
              rhoval = rhoval + rho(n,R_A)
              gamma  = gamma + rgamma(n,G_AA)
            endif
            if (rho(n,R_B).ge.0.5d0*tol_rho) then
              rhoval = rhoval + rho(n,R_B)
              gamma  = gamma + rgamma(n,G_BB)
              if (rho(n,R_A).ge.0.5d0*tol_rho) then
                gamma  = gamma + 2.0d0*rgamma(n,G_AB)
              endif
            endif
            if (rhoval.lt.tol_rho) goto 20
            arho=1.d0/rhoval
            rho13  = abs(rhoval)**ONE3
            rho43  = rhoval*rho13
            rho76  = abs(rhoval)**SEV6
            rs = rsfact/rho13
            rs2 = rs*rs
            rs3 = rs2*rs
c           d1rs = -ONE3*rs*arho
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
c           d2rs = -FOUR3*d1rs*arho
#endif
#if defined(THIRD_DERIV)
c           d3rs = -SEVEN3*d2rs*arho 
#endif
c           gamma = rgamma(n,G_AA)+rgamma(n,G_BB)+2.0d0*rgamma(n,G_AB)
c           gamma = delrho(n,1,1)*delrho(n,1,1) +
c    &              delrho(n,2,1)*delrho(n,2,1) +
c    &              delrho(n,3,1)*delrho(n,3,1) +
c    &              delrho(n,1,2)*delrho(n,1,2) +
c    &              delrho(n,2,2)*delrho(n,2,2) +
c    &              delrho(n,3,2)*delrho(n,3,2) +
c    &        2.d0*(delrho(n,1,1)*delrho(n,1,2) +
c    &              delrho(n,2,1)*delrho(n,2,2) +
c    &              delrho(n,3,1)*delrho(n,3,2))
            if (gamma.gt.tol_rho*tol_rho) then
              gam12 = sqrt(gamma)
            else
              gam12 = 0.0d0
            endif
            zeta = (rho(n,R_A) - rho(n,R_B))*arho
            if(zeta.le.-1d0) zeta=-1d0
            if(zeta.ge.1d0) zeta=1d0
c           d1z(1) =  (1.d0 - zeta)*arho
c           d1z(2) = -(1.d0 + zeta)*arho
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
c           rrho2 = 2.d0*arho*arho
c           1 = aa, 2 = ab, 3 = bb
c           d2z(1) =-rrho2*(1.d0-zeta)
c           d2z(2) = rrho2*zeta
c           d2z(3) = rrho2*(1.d0+zeta)
#endif
#if defined(THIRD_DERIV)
c           d3rs = -SEVEN3*d2rs*arho
c           if ((1.d0-zeta).lt.tol_rho) then
c             d3fz = (1.d0+zeta)**(-FIVE3)
c           else if ((1.d0+zeta).lt.tol_rho) then
c             d3fz = (1.d0-zeta)**(-FIVE3)
c           else
c             d3fz = (1.d0+zeta)**(-FIVE3) + (1.d0-zeta)**(-FIVE3)
c           end if
c           d3fz = -d3fz*TWO3*ONE3*FOUR3/(2.d0**FOUR3-2.d0)
c
c           rrho3 = rrho2*arho
c
c 1 = aaa, 2 = aab, 3 = abb, 4 = bbb
c           d3z(1) = 3.0d0*rrho3*(1.0d0 - zeta)
c           d3z(2) = rrho3*(1.0d0 - 3.0d0*zeta)
c           d3z(3) = -rrho3*(1.0d0 + 3.0d0*zeta)
c           d3z(4) = -3.0d0*rrho3*(1.0d0 + zeta)
#endif
c
c           d(zeta)
c
            dt12 = 0.0d0
            if (ONE+zeta.gt.1.0d-10) then
              dt12 = dt12 + (0.5d0*(ONE+zeta))**FIVE3
            endif
            if (ONE-zeta.gt.1.0d-10) then
              dt12 = dt12 + (0.5d0*(ONE-zeta))**FIVE3
            endif
c           d1dt12 = FIVE3*0.5d0*(
c    &           ((ONE+zeta)*.5d0)**TWO3 - ((ONE-zeta)*.5d0)**TWO3 )
            d = 2.d0**ONE3*sqrt(dt12)
            dm1 = 1.d0/d
c           adp = 0.5d0*d/dt12*d1dt12
c           d1d(1) = adp*d1z(1)
c           d1d(2) = adp*d1z(2)
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
c           if ((1.d0-zeta).lt.tol_rho) then
c             d2dt12 = FIVE3*TWO3*0.25d0*(((ONE+zeta)*.5d0)**(-ONE3))
c           else if ((1.d0+zeta).lt.tol_rho) then
c             d2dt12 = FIVE3*TWO3*0.25d0*(((ONE-zeta)*.5d0)**(-ONE3))
c           else
c             d2dt12 = FIVE3*TWO3*0.25d0*(
c    &         ((ONE+zeta)*.5d0)**(-ONE3) + ((ONE-zeta)*.5d0)**(-ONE3) )
c           end if
c
c           dpp =-0.5d0*adp/dt12*d1dt12
c    &        + 2.d0**(-TWO3)*d2dt12/dsqrt(dt12)
c           d2d(1) = dpp*d1z(1)*d1z(1) + adp*d2z(1)
c           d2d(2) = dpp*d1z(1)*d1z(2) + adp*d2z(2)
c           d2d(3) = dpp*d1z(2)*d1z(2) + adp*d2z(3)
#endif
#if defined(THIRD_DERIV)
c           call errquit("nwxc_c_perdew86: no 3rd derivatives",0,0)
#endif
c         
c           C(n)
c         
            anum = fff+alpha*rs+beta*rs2
            aden = 1.d0+pgamma*rs+delta*rs2+beta10*rs3
            Cn = zzz + anum/aden
c           d1anum = alpha + 2d0*beta*rs
c           d1aden = pgamma + 2d0*delta*rs + 3d0*beta10*rs2
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
c           d2anum = 2d0*beta
c           d2aden = 2d0*delta + 6d0*beta10*rs
#endif
c     First compute rs derivative
c           d1Cn = d1anum/aden - anum*d1aden/aden**2
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
c           d2Cn = d2anum/aden - (2d0*d1anum*d1aden+anum*d2aden)/aden**2
c    &           + 2d0*anum*d1aden**2/aden**3
#endif
c     Convert to rho derivative
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
c           d2Cn = d2Cn*d1rs*d1rs + d1Cn*d2rs
#endif
c           d1Cn = d1Cn*d1rs
c         
c           phi(n,gradn)
c         
            expfac = 0.d0
            phi = (pfff*ftilde)*(CINF/Cn)*gam12/rho76
            if (phi.lt.EXPTOL) expfac = exp(-phi)
c           dlnphi = -(d1Cn/Cn + SEV6/rhoval)
c           d1phi(1) = phi*dlnphi
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
c           d2phi(1) = d1phi(1)*dlnphi
c    &               + phi*((d1Cn/Cn)**2 - d2Cn/Cn + SEV6/rhoval**2)
#endif
c         
c           functional
c         
            func = expfac*Cn*gamma/rho43*dm1
c           t = d1Cn/Cn - (d1phi(1) + FOUR3/rhoval)
c           dlnfrho(1) = t - dm1*d1d(1)
c           dlnfrho(2) = t - dm1*d1d(2)
c           d1f(1) = dlnfrho(1)*func
c           d1f(2) = dlnfrho(2)*func
c           Amat(n,D1_RA) = Amat(n,D1_RA) + d1f(1)*wght
c           Amat(n,D1_RB) = Amat(n,D1_RB) + d1f(2)*wght
c           if (gam12.gt.TOLL)then
c              d1phi(2) = phi / (2d0*gamma)
c              dlnfgam = 1d0/gamma - d1phi(2)
c              d1f(3) = func*dlnfgam
c              Cmat(n,D1_GAA) = Cmat(n,D1_GAA) + d1f(3)*wght
c              Cmat(n,D1_GAB) = Cmat(n,D1_GAB) + d1f(3)*2D0*wght
c              Cmat(n,D1_GBB) = Cmat(n,D1_GBB) + d1f(3)*wght
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
c              d2phi(2) = d1phi(2)*dlnphi
c              d2phi(3) =-d1phi(2)/(2d0*gamma)
c
c              t = -d2phi(1) - (d1Cn/Cn)**2 + d2Cn/Cn + FOUR3/rhoval**2
c              Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA)
c    &              + (d1f(1)*dlnfrho(1)
c    &              + func*(d1d(1)*d1d(1)*dm1**2-d2d(1)*dm1+t))*wght
c              Amat2(n,D2_RA_RB) = Amat2(n,D2_RA_RB)
c    &              + (d1f(1)*dlnfrho(2)
c    &              + func*(d1d(1)*d1d(2)*dm1**2-d2d(2)*dm1+t))*wght
c              Amat2(n,D2_RB_RB) = Amat2(n,D2_RB_RB)
c    &              + (d1f(2)*dlnfrho(2)
c    &              + func*(d1d(2)*d1d(2)*dm1**2-d2d(3)*dm1+t))*wght
c
c              t = (d1f(1)*dlnfgam - func*d2phi(2))*wght
c              Cmat2(n,D2_RA_GAA) = Cmat2(n,D2_RA_GAA) + t
c              Cmat2(n,D2_RA_GAB) = Cmat2(n,D2_RA_GAB) + t*2d0
c              Cmat2(n,D2_RA_GBB) = Cmat2(n,D2_RA_GBB) + t
c              t = (d1f(2)*dlnfgam - func*d2phi(2))*wght
c              Cmat2(n,D2_RB_GAA) = Cmat2(n,D2_RB_GAA) + t
c              Cmat2(n,D2_RB_GAB) = Cmat2(n,D2_RB_GAB) + t*2d0
c              Cmat2(n,D2_RB_GBB) = Cmat2(n,D2_RB_GBB) + t
c
c              t = (d1f(3)*dlnfgam - func*(1d0/gamma**2+d2phi(3)))*wght
c              Cmat2(n,D2_GAA_GAA) = Cmat2(n,D2_GAA_GAA) + t
c              Cmat2(n,D2_GAA_GBB) = Cmat2(n,D2_GAA_GBB) + t
c              Cmat2(n,D2_GBB_GBB) = Cmat2(n,D2_GBB_GBB) + t
c              Cmat2(n,D2_GAA_GAB) = Cmat2(n,D2_GAA_GAB) + t*2d0
c              Cmat2(n,D2_GAB_GBB) = Cmat2(n,D2_GAB_GBB) + t*2d0
c              Cmat2(n,D2_GAB_GAB) = Cmat2(n,D2_GAB_GAB) + t*4d0
#endif
c           endif
            ffunc(n)=ffunc(n)+func*wght
   20    continue
      endif
      return
      end
#ifndef NWAD_PRINT
#define NWAD_PRINT
c
c     Compile source again for Maxima
c
#include "nwxc_c_perdew86.F"
#endif
#ifndef SECOND_DERIV
#define SECOND_DERIV
c
c     Compile source again for the 2nd derivative case
c
#include "nwxc_c_perdew86.F"
#endif
#ifndef THIRD_DERIV
#define THIRD_DERIV
c
c     Compile source again for the 3rd derivative case
c
#include "nwxc_c_perdew86.F"
#endif
#undef NWAD_PRINT
C>
C> @}
