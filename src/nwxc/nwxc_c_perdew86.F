#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
C> \ingroup nwxc
C> @{
C>
C> \file nwxc_c_perdew86.F
C> The Perdew correlation functional of 1986
C>
C> @}
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
C>     <a href="http://dx.doi.org/10.1103/PhysRevB.33.8822">
C>     10.1103/PhysRevB.33.8822</a>.
C>
C> [2] P. Mlynarski, D.R. Salahub,
C>     "Self-consistent implementation of nonlocal exchange and
C>     correlation in a Gaussian density-functional method",
C>     Phys. Rev. B <b>43</b>, 1399–1410 (1991), DOI:
C>     <a href="http://dx.doi.org/10.1103/PhysRevB.43.1399">
C>     10.1103/PhysRevB.43.1399</a>.
C>
C> [3] J.P. Perdew,
C>     "Erratum: Density-functional approximation for the correlation
C>     energy of the inhomogeneous electron gas", Phys. Rev. B
C>     <b>34</b>, 7406–7406 (1986), DOI:
C>     <a href="http://dx.doi.org/10.1103/PhysRevB.34.7406">
C>     10.1103/PhysRevB.34.7406</a>.
C>
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
      Subroutine nwxc_c_perdew86(tol_rho, ipol, nq, wght, rho, rgamma, 
     &                           ffunc, Amat, Cmat)
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
      Subroutine nwxc_c_perdew86_d2(tol_rho, ipol, nq, wght,
     &                              rho, rgamma, ffunc, 
     &                              Amat, Amat2, Cmat, Cmat2)
#else
      Subroutine nwxc_c_perdew86_d3(tol_rho, ipol, nq, wght,
     &                              rho, rgamma, ffunc, 
     &                              Amat, Amat2, Amat3,
     &                              Cmat, Cmat2, Cmat3)
#endif
c
c$Id$
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
      double precision rho(nq,*)    !< [Input] The density
c
c     Charge Density Gradient
c
      double precision rgamma(nq,*) !< [Input] The norm of the density gradients
c
c     Sampling Matrices for the XC Potential
c
      double precision ffunc(nq)    !< [Output] The value of the functional
      double precision Amat(nq,*)   !< [Output] The derivative wrt rho
      double precision Cmat(nq,*)   !< [Output] The derivative wrt rgamma
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
c
c     Sampling Matrices for the XC Kernel
c
      double precision Amat2(nq,*)  !< [Output] The 2nd derivative wrt rho
      double precision Cmat2(nq,*)  !< [Output] The 2nd derivative wrt rgamma
                                    !< and possibly rho
#endif
#if defined(THIRD_DERIV)
c
c     Sampling Matrices for the XC Kernel
c
      double precision Amat3(nq,*)  !< [Output] The 3rd derivative wrt rho
      double precision Cmat3(nq,*)  !< [Output] The 3rd derivative wrt rgamma
                                    !< and possibly rho
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
      double precision rsfact, rs, rs2, rs3
      double precision rhoval, rho13, rho43, rho76, arho
      double precision d1rs
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
      double precision d2rs
#endif
#if defined(THIRD_DERIV)
      double precision d3rs
#endif
      double precision gamma, gam12
      double precision anum, aden, d1anum, d1aden, Cn, d1Cn,
     &     expfac, phi, d1phi(2), dlnphi, func, d1f(3),
     &     dlnfrho(2), dlnfgam
      double precision zeta, d1z(2), d, dm1, adp, d1d(2), t,
     &     dt12, d1dt12
      double precision aden2
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
      double precision d2anum, d2aden, rrho2, d2z(3), dpp, d2d(3),
     &     d2phi(3), d2dt12, d2Cn
      double precision aden3
      double precision arho2
      double precision d2lnphi
      double precision d2f(3)
      double precision d2lnfrho(3), d2lnfrg(2), d2lnfgam
#endif
#if defined(THIRD_DERIV)
      double precision d3lnphi
      double precision d3anum, d3aden, d3Cn, d3phi(4)
      double precision d3lnfrho(4), d3lnfgam
      double precision d3f(3)
      double precision aden4
      double precision arho3
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
            rho13 = abs(rhoval)**ONE3
            rho43 = rhoval*rho13
            rho76 = abs(rhoval)**SEV6
            rs = rsfact/rho13
            rs2 = rs*rs
            rs3 = rs2*rs
            d1rs = -ONE3*rs*arho
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
            d2rs = -FOUR3*d1rs*arho
#endif
#if defined(THIRD_DERIV)
            d3rs = -SEVEN3*d2rs*arho 
#endif
            gamma = rgamma(n,G_TT)
c           gamma = delrho(n,1,1)*delrho(n,1,1) +
c    &              delrho(n,2,1)*delrho(n,2,1) +
c    &              delrho(n,3,1)*delrho(n,3,1)
            gam12 = sqrt(gamma)
c         
c           C(n)
c         
            anum = fff+alpha*rs+beta*rs2
            aden = 1.d0+pgamma*rs+delta*rs2+beta10*rs3
            Cn = zzz + anum/aden
            d1anum = alpha + 2d0*beta*rs
            d1aden = pgamma + 2d0*delta*rs + 3d0*beta10*rs2
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
            d2anum = 2d0*beta
            d2aden = 2d0*delta + 6d0*beta10*rs
#endif
#if defined(THIRD_DERIV)
            d3anum = 0.0d0
            d3aden = 6.0d0*beta10 
#endif
c     First compute rs derivative
            aden2 = aden*aden
            d1Cn = d1anum/aden - anum*d1aden/aden2
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
            aden3 = aden2*aden
            d2Cn = d2anum/aden - (2d0*d1anum*d1aden+anum*d2aden)/aden2
     &           + 2d0*anum*d1aden**2/aden3
#endif
#if defined(THIRD_DERIV)
            aden4 = aden3*aden
c
            d3Cn = -( 3.0d0*d2anum*d1aden + 3.0d0*d1anum*d2aden
     1              + anum*d3aden )/aden2
     2           + 6.0d0*( d1anum*d1aden**2
     3                   + anum*d2aden*d1aden )/aden3
     4           - 6.0d0*anum*d1aden**3/aden4
#endif
c     Convert to rho derivative
#if defined(THIRD_DERIV)
            d3Cn = d3Cn*d1rs*d1rs*d1rs
     1           + 3.0d0*d2Cn*d2rs*d1rs
     2           + d1Cn*d3rs
#endif
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
            d2Cn = d2Cn*d1rs*d1rs + d1Cn*d2rs
#endif
            d1Cn = d1Cn*d1rs
c         
c           phi(n,gradn)
c         
            expfac = 0.d0
            phi = (pfff*ftilde)*(CINF/Cn)*gam12/rho76
            if (phi.lt.EXPTOL) expfac = exp(-phi)
            dlnphi = -(d1Cn/Cn + SEV6/rhoval)
            d1phi(1) = phi*dlnphi
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
            arho2 = arho*arho
            d2lnphi = (d1Cn/Cn)**2 - d2Cn/Cn + SEV6*arho2
            d2phi(1) = d1phi(1)*dlnphi + phi*d2lnphi
c           d2phi(1) = d1phi(1)*dlnphi
c    &               + phi*((d1Cn/Cn)**2 - d2Cn/Cn + SEV6/rhoval**2)
#endif
#if defined(THIRD_DERIV)
            arho3 = arho2*arho
c
            d3lnphi = -2.0d0*(d1Cn/Cn)**3
     1              + 3.0d0*(d2Cn/Cn)*(d1Cn/Cn)
     2              - d3Cn/Cn
     3              - SEVEN3*arho3
            d3phi(1) = d2phi(1)*dlnphi
     1               + 2.0d0*d1phi(1)*d2lnphi
     2               + phi*d3lnphi
#endif
c         
c           functional
c         
            func = expfac*Cn*gamma/rho43
            dlnfrho(1) = d1Cn/Cn - (d1phi(1) + FOUR3/rhoval)
            d1f(1) = dlnfrho(1)*func
            Amat(n,D1_RA) = Amat(n,D1_RA) + d1f(1)*wght
            if (gam12.gt.TOLL)then
               d1phi(2) = phi / (2d0*gamma)
               dlnfgam = 1d0/gamma - d1phi(2)
               d1f(3) = func*dlnfgam
               Cmat(n,D1_GAA) = Cmat(n,D1_GAA) + d1f(3)*wght
               Cmat(n,D1_GAB) = Cmat(n,D1_GAB) + d1f(3)*2D0*wght
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
               d2phi(2) = d1phi(2)*dlnphi
               d2phi(3) =-d1phi(2)/(2d0*gamma)
c!!! Which of the following are actually needed for restricted?
c!!! Should treat derivatives of d as zero? d is a constant?
c Daniel (11-19-12): d is a constant (it equals 1) for a restricted 
c calculation, since there is no spin-polarization.  Thus, the
c derivatives are zero.
               d2lnfrho(1) = -d2phi(1) - (d1Cn/Cn)**2 + d2Cn/Cn
     1                     + FOUR3*arho2
c
               d2f(1) = d1f(1)*dlnfrho(1)
     1                + func*d2lnfrho(1)
c
               t = d2f(1)*wght
c
               Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA) + t
               Amat2(n,D2_RA_RB) = Amat2(n,D2_RA_RB) + t
     &              + (d1f(1)*dlnfrho(1)
     &              + func*t)*wght
#if 0
               Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA)
     &              + (d1f(1)*dlnfrho(1)
     &              + func*(d1d(1)*d1d(1)*dm1**2-d2d(1)*dm1+t))*wght
               Amat2(n,D2_RA_RB) = Amat2(n,D2_RA_RB)
     &              + (d1f(1)*dlnfrho(2)
     &              + func*(d1d(1)*d1d(2)*dm1**2-d2d(2)*dm1+t))*wght
#endif
c rg terms
               d2lnfrg(1) = -d2phi(2)
               d2f(2) = (d1f(1)*dlnfgam + func*d2lnfrg(1))
               t = d2f(2)*wght
c
               Cmat2(n,D2_RA_GAA) = Cmat2(n,D2_RA_GAA) + t
               Cmat2(n,D2_RA_GAB) = Cmat2(n,D2_RA_GAB) + t*2d0
               Cmat2(n,D2_RA_GBB) = Cmat2(n,D2_RA_GBB) + t
c gg terms
               d2lnfgam = -1.0d0/gamma**2 - d2phi(3)
               d2f(3) = d1f(3)*dlnfgam + func*d2lnfgam
               t = d2f(3)*wght
c
               Cmat2(n,D2_GAA_GAA) = Cmat2(n,D2_GAA_GAA) + t
               Cmat2(n,D2_GAA_GBB) = Cmat2(n,D2_GAA_GBB) + t
               Cmat2(n,D2_GAA_GAB) = Cmat2(n,D2_GAA_GAB) + t*2d0
               Cmat2(n,D2_GAB_GAB) = Cmat2(n,D2_GAB_GAB) + t*4d0
#endif
#if defined(THIRD_DERIV)
c rrr terms
               d3lnfrho(1) = -d3phi(1)
     1                     + 2.0d0*(d1Cn/Cn)**3
     2                     - 3.0d0*(d2Cn/Cn)*(d1Cn/Cn)
     3                     + d3Cn/Cn
     4                     - EIGHT3*arho3
c
               d3f(1) = d2f(1)*dlnfrho(1)
     1                + 2.0d0*d1f(1)*d2lnfrho(1)
     2                + func*d3lnfrho(1)
c
               t = d3f(1)*wght
c
               Amat3(n,D3_RA_RA_RA) = Amat3(n,D3_RA_RA_RA) + t
               Amat3(n,D3_RA_RA_RB) = Amat3(n,D3_RA_RA_RB) + t
               Amat3(n,D3_RA_RB_RB) = Amat3(n,D3_RA_RB_RB) + t
c rrg terms
               d3phi(2) = d2phi(2)*dlnphi + d1phi(2)*d2lnphi
c
               t = ( d2f(2)*dlnfrho(1)
     1             - d1f(1)*d2phi(2)
     2             + d1f(3)*d2lnfrho(1)
     3             - func*d3phi(2) )*wght
c
               Cmat3(n,D3_RA_RA_GAA) = Cmat3(n,D3_RA_RA_GAA) + t
               Cmat3(n,D3_RA_RA_GAB) = Cmat3(n,D3_RA_RA_GAB) + t*2.0d0
               Cmat3(n,D3_RA_RA_GBB) = Cmat3(n,D3_RA_RA_GBB) + t
               Cmat3(n,D3_RA_RB_GAA) = Cmat3(n,D3_RA_RB_GAA) + t
               Cmat3(n,D3_RA_RB_GAB) = Cmat3(n,D3_RA_RB_GAB) + t*2.0d0
               Cmat3(n,D3_RA_RB_GBB) = Cmat3(n,D3_RA_RB_GBB) + t
c rgg terms
               d3phi(3) = -d2phi(3)*dlnphi
c
               t = ( d2f(2)*dlnfgam
     1             + d1f(1)*d2lnfgam
     2             + d1f(3)*d2lnfrg(1)
     3             + func*d3phi(3) )*wght
c
               Cmat3(n,D3_RA_GAA_GAA) = Cmat3(n,D3_RA_GAA_GAA) + t
               Cmat3(n,D3_RA_GAA_GAB) = Cmat3(n,D3_RA_GAA_GAB) + t*2.0d0
               Cmat3(n,D3_RA_GAA_GBB) = Cmat3(n,D3_RA_GAA_GBB) + t
               Cmat3(n,D3_RA_GAB_GAB) = Cmat3(n,D3_RA_GAB_GAB) + t*4.0d0
               Cmat3(n,D3_RA_GAB_GBB) = Cmat3(n,D3_RA_GAB_GBB) + t*2.0d0
               Cmat3(n,D3_RA_GBB_GBB) = Cmat3(n,D3_RA_GBB_GBB) + t
c ggg terms              
               d3phi(4) = -3.0d0*d2phi(3)/(2.0d0*gamma)
               d3lnfgam = 2.0d0/gamma**3 - d3phi(4)
c
               t = ( d2f(3)*dlnfgam
     1             + 2.0d0*d1f(3)*d2lnfgam
     2             + func*d3lnfgam )*wght
c
               Cmat3(n,D3_GAA_GAA_GAA) = Cmat3(n,D3_GAA_GAA_GAA) + t
               Cmat3(n,D3_GAA_GAA_GAB) = Cmat3(n,D3_GAA_GAA_GAB) 
     1                                 + t*2.0d0
               Cmat3(n,D3_GAA_GAA_GBB) = Cmat3(n,D3_GAA_GAA_GBB) + t
               Cmat3(n,D3_GAA_GAB_GAB) = Cmat3(n,D3_GAA_GAB_GAB) 
     1                                 + t*4.0d0
               Cmat3(n,D3_GAA_GAB_GBB) = Cmat3(n,D3_GAA_GAB_GBB) 
     1                                 + t*2.0d0
               Cmat3(n,D3_GAA_GBB_GBB) = Cmat3(n,D3_GAA_GBB_GBB) + t
               Cmat3(n,D3_GAB_GAB_GAB) = Cmat3(n,D3_GAB_GAB_GAB) 
     1                                 + t*8.0d0
#endif
            endif
            ffunc(n)=ffunc(n)+func*wght
   10    continue
      else
c
c        ======> SPIN-UNRESTRICTED <======
c
         do 20 n = 1, nq
            rhoval = rho(n,R_A)+rho(n,R_B)
            if (rhoval.lt.tol_rho) goto 20
            arho=1.d0/rhoval
            rho13  = abs(rhoval)**ONE3
            rho43  = rhoval*rho13
            rho76  = abs(rhoval)**SEV6
            rs = rsfact/rho13
            rs2 = rs*rs
            rs3 = rs2*rs
            d1rs = -ONE3*rs*arho
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
            d2rs = -FOUR3*d1rs*arho
#endif
#if defined(THIRD_DERIV)
            d3rs = -SEVEN3*d2rs*arho 
#endif
            gamma = rgamma(n,G_AA)+rgamma(n,G_BB)+2.0d0*rgamma(n,G_AB)
c           gamma = delrho(n,1,1)*delrho(n,1,1) +
c    &              delrho(n,2,1)*delrho(n,2,1) +
c    &              delrho(n,3,1)*delrho(n,3,1) +
c    &              delrho(n,1,2)*delrho(n,1,2) +
c    &              delrho(n,2,2)*delrho(n,2,2) +
c    &              delrho(n,3,2)*delrho(n,3,2) +
c    &        2.d0*(delrho(n,1,1)*delrho(n,1,2) +
c    &              delrho(n,2,1)*delrho(n,2,2) +
c    &              delrho(n,3,1)*delrho(n,3,2))
            gam12 = sqrt(gamma)
            zeta = (rho(n,R_A) - rho(n,R_B))*arho
            if(zeta.lt.-1d0) zeta=-1d0
            if(zeta.gt.1d0) zeta=1d0
            d1z(1) =  (1.d0 - zeta)*arho
            d1z(2) = -(1.d0 + zeta)*arho
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
            rrho2 = 2.d0*arho*arho
c           1 = aa, 2 = ab, 3 = bb
            d2z(1) =-rrho2*(1.d0-zeta)
            d2z(2) = rrho2*zeta
            d2z(3) = rrho2*(1.d0+zeta)
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
            dt12 = ((ONE+zeta)*.5d0)**FIVE3 + ((ONE-zeta)*.5d0)**FIVE3
            d1dt12 = FIVE3*0.5d0*(
     &           ((ONE+zeta)*.5d0)**TWO3 - ((ONE-zeta)*.5d0)**TWO3 )
            d = 2.d0**ONE3*dsqrt(dt12)
            dm1 = 1.d0/d
            adp = 0.5d0*d/dt12*d1dt12
            d1d(1) = adp*d1z(1)
            d1d(2) = adp*d1z(2)
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
            if ((1.d0-zeta).lt.tol_rho) then
              d2dt12 = FIVE3*TWO3*0.25d0*(((ONE+zeta)*.5d0)**(-ONE3))
            else if ((1.d0+zeta).lt.tol_rho) then
              d2dt12 = FIVE3*TWO3*0.25d0*(((ONE-zeta)*.5d0)**(-ONE3))
            else
              d2dt12 = FIVE3*TWO3*0.25d0*(
     &         ((ONE+zeta)*.5d0)**(-ONE3) + ((ONE-zeta)*.5d0)**(-ONE3) )
            end if
c
            dpp =-0.5d0*adp/dt12*d1dt12
     &        + 2.d0**(-TWO3)*d2dt12/dsqrt(dt12)
            d2d(1) = dpp*d1z(1)*d1z(1) + adp*d2z(1)
            d2d(2) = dpp*d1z(1)*d1z(2) + adp*d2z(2)
            d2d(3) = dpp*d1z(2)*d1z(2) + adp*d2z(3)
#endif
#if defined(THIRD_DERIV)
            call errquit("nwxc_c_perdew86: no 3rd derivatives",0,0)
#endif
c         
c           C(n)
c         
            anum = fff+alpha*rs+beta*rs2
            aden = 1.d0+pgamma*rs+delta*rs2+beta10*rs3
            Cn = zzz + anum/aden
            d1anum = alpha + 2d0*beta*rs
            d1aden = pgamma + 2d0*delta*rs + 3d0*beta10*rs2
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
            d2anum = 2d0*beta
            d2aden = 2d0*delta + 6d0*beta10*rs
#endif
c     First compute rs derivative
            d1Cn = d1anum/aden - anum*d1aden/aden**2
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
            d2Cn = d2anum/aden - (2d0*d1anum*d1aden+anum*d2aden)/aden**2
     &           + 2d0*anum*d1aden**2/aden**3
#endif
c     Convert to rho derivative
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
            d2Cn = d2Cn*d1rs*d1rs + d1Cn*d2rs
#endif
            d1Cn = d1Cn*d1rs
c         
c           phi(n,gradn)
c         
            expfac = 0.d0
            phi = (pfff*ftilde)*(CINF/Cn)*gam12/rho76
            if (phi.lt.EXPTOL) expfac = exp(-phi)
            dlnphi = -(d1Cn/Cn + SEV6/rhoval)
            d1phi(1) = phi*dlnphi
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
            d2phi(1) = d1phi(1)*dlnphi
     &               + phi*((d1Cn/Cn)**2 - d2Cn/Cn + SEV6/rhoval**2)
#endif
c         
c           functional
c         
            func = expfac*Cn*gamma/rho43*dm1
            t = d1Cn/Cn - (d1phi(1) + FOUR3/rhoval)
            dlnfrho(1) = t - dm1*d1d(1)
            dlnfrho(2) = t - dm1*d1d(2)
            d1f(1) = dlnfrho(1)*func
            d1f(2) = dlnfrho(2)*func
            Amat(n,D1_RA) = Amat(n,D1_RA) + d1f(1)*wght
            Amat(n,D1_RB) = Amat(n,D1_RB) + d1f(2)*wght
            if (gam12.gt.TOLL)then
               d1phi(2) = phi / (2d0*gamma)
               dlnfgam = 1d0/gamma - d1phi(2)
               d1f(3) = func*dlnfgam
               Cmat(n,D1_GAA) = Cmat(n,D1_GAA) + d1f(3)*wght
               Cmat(n,D1_GAB) = Cmat(n,D1_GAB) + d1f(3)*2D0*wght
               Cmat(n,D1_GBB) = Cmat(n,D1_GBB) + d1f(3)*wght
#if defined(SECOND_DERIV) || defined(THIRD_DERIV)
               d2phi(2) = d1phi(2)*dlnphi
               d2phi(3) =-d1phi(2)/(2d0*gamma)
c
               t = -d2phi(1) - (d1Cn/Cn)**2 + d2Cn/Cn + FOUR3/rhoval**2
               Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA)
     &              + (d1f(1)*dlnfrho(1)
     &              + func*(d1d(1)*d1d(1)*dm1**2-d2d(1)*dm1+t))*wght
               Amat2(n,D2_RA_RB) = Amat2(n,D2_RA_RB)
     &              + (d1f(1)*dlnfrho(2)
     &              + func*(d1d(1)*d1d(2)*dm1**2-d2d(2)*dm1+t))*wght
               Amat2(n,D2_RB_RB) = Amat2(n,D2_RB_RB)
     &              + (d1f(2)*dlnfrho(2)
     &              + func*(d1d(2)*d1d(2)*dm1**2-d2d(3)*dm1+t))*wght
c
               t = (d1f(1)*dlnfgam - func*d2phi(2))*wght
               Cmat2(n,D2_RA_GAA) = Cmat2(n,D2_RA_GAA) + t
               Cmat2(n,D2_RA_GAB) = Cmat2(n,D2_RA_GAB) + t*2d0
               Cmat2(n,D2_RA_GBB) = Cmat2(n,D2_RA_GBB) + t
               t = (d1f(2)*dlnfgam - func*d2phi(2))*wght
               Cmat2(n,D2_RB_GAA) = Cmat2(n,D2_RB_GAA) + t
               Cmat2(n,D2_RB_GAB) = Cmat2(n,D2_RB_GAB) + t*2d0
               Cmat2(n,D2_RB_GBB) = Cmat2(n,D2_RB_GBB) + t
c
               t = (d1f(3)*dlnfgam - func*(1d0/gamma**2+d2phi(3)))*wght
               Cmat2(n,D2_GAA_GAA) = Cmat2(n,D2_GAA_GAA) + t
               Cmat2(n,D2_GAA_GBB) = Cmat2(n,D2_GAA_GBB) + t
               Cmat2(n,D2_GBB_GBB) = Cmat2(n,D2_GBB_GBB) + t
               Cmat2(n,D2_GAA_GAB) = Cmat2(n,D2_GAA_GAB) + t*2d0
               Cmat2(n,D2_GAB_GBB) = Cmat2(n,D2_GAB_GBB) + t*2d0
               Cmat2(n,D2_GAB_GAB) = Cmat2(n,D2_GAB_GAB) + t*4d0
#endif
            endif
            ffunc(n)=ffunc(n)+func*wght
   20    continue
      endif
      return
      end
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
C>
C> @}

