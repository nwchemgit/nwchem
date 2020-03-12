#ifndef SECOND_DERIV
C> \ingroup nwpwxc
C> @{
C>
C> \file nwpwxc_c_perdew86.F
C> The Perdew correlation functional of 1986
C>
C> @}
#endif
C> \ingroup nwpwxc_priv
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
#ifndef SECOND_DERIV
      Subroutine nwpwxc_c_perdew86(tol_rho, ipol, nq, wght, rho, rgamma, 
     &                           ffunc, Amat, Cmat)
#else
      Subroutine nwpwxc_c_perdew86_d2(tol_rho, ipol, nq, wght,
     &                              rho, rgamma, ffunc, 
     &                              Amat, Amat2, Cmat, Cmat2)
#endif
c
c$Id: xc_perdew86.F 22803 2012-08-30 05:26:01Z niri $
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
      double precision ffunc(nq)    !< [Output] The value of the functional
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
      double precision TOLL, EXPTOL, alpha, beta, pgamma, delta, 
     &                 beta10, ftilde, zzz, fff, pfff, CINF, ONE, 
     &                 ONE3, THREE, FOUR3, SEV6, FIVE3, 
     &                 TWO3, FIVE6, pi
      Parameter (TOLL = 1.D-40, EXPTOL = 80.d0)
      Parameter (alpha = 0.023266D0, beta  =  7.389D-6, 
     &   pgamma = 8.723d0, delta = 0.472d0,  beta10 = 10000.d0*beta)
      parameter (ftilde = 0.11d0, zzz = 0.001667d0, fff = 0.002568d0)
      parameter(pfff = 1.745d0, CINF = zzz+fff)
      Parameter (ONE = 1.D0, ONE3 = 1.d0/3.d0, THREE = 3.d0)
      Parameter (FOUR3 = 4.D0/3.D0, SEV6 = 7.d0/6.d0)
      parameter (FIVE3 = 5.d0/3.d0, TWO3 = 2.d0/3.d0, FIVE6 = 5.d0/6.d0)
c     parameter (pi = 3.1415926535897932385d0)
c     
c     Mlynarski Salahub PRB 43, 1399 (1991)
c     
      integer n
      double precision rsfact, rs, rs2, rs3
      double precision rhoval, rho13, rho43, rho76, arho
      double precision d1rs
#ifdef SECOND_DERIV
      double precision d2rs
#endif
      double precision gamma, gam12
      double precision anum, aden, d1anum, d1aden, Cn, d1Cn,
     &     expfac, phi, d1phi(2), dlnphi, func, d1f(3),
     &     dlnfrho(2), dlnfgam
      double precision zeta, d1z(2), d, dm1, adp, d1d(2), t,
     &     dt12, d1dt12
#ifdef SECOND_DERIV
      double precision d2anum, d2aden, rrho2, d2z(3), dpp, d2d(3),
     &     d2phi(3), d2dt12, d2Cn
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
#ifdef SECOND_DERIV
            d2rs = -FOUR3*d1rs*arho
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
#ifdef SECOND_DERIV
            d2anum = 2d0*beta
            d2aden = 2d0*delta + 6d0*beta10*rs
#endif
c     First compute rs derivative
            d1Cn = d1anum/aden - anum*d1aden/aden**2
#ifdef SECOND_DERIV
            d2Cn = d2anum/aden - (2d0*d1anum*d1aden+anum*d2aden)/aden**2
     &           + 2d0*anum*d1aden**2/aden**3
#endif
c     Convert to rho derivative
#ifdef SECOND_DERIV
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
#ifdef SECOND_DERIV
            d2phi(1) = d1phi(1)*dlnphi
     &               + phi*((d1Cn/Cn)**2 - d2Cn/Cn + SEV6/rhoval**2)
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
#ifdef SECOND_DERIV
               d2phi(2) = d1phi(2)*dlnphi
               d2phi(3) =-d1phi(2)/(2d0*gamma)
c!!! Which of the following are actually needed for restricted?
c!!! Should treat derivatives of d as zero? d is a constant?
               t = -d2phi(1) - (d1Cn/Cn)**2 + d2Cn/Cn + FOUR3/rhoval**2
               Amat2(n,D2_RA_RA) = Amat2(n,D2_RA_RA)
     &              + (d1f(1)*dlnfrho(1)
     &              + func*t)*wght
               Amat2(n,D2_RA_RB) = Amat2(n,D2_RA_RB)
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
c
               t = (d1f(1)*dlnfgam - func*d2phi(2))*wght
               Cmat2(n,D2_RA_GAA) = Cmat2(n,D2_RA_GAA) + t
               Cmat2(n,D2_RA_GAB) = Cmat2(n,D2_RA_GAB) + t*2d0
               Cmat2(n,D2_RA_GBB) = Cmat2(n,D2_RA_GBB) + t
c
               t = (d1f(3)*dlnfgam - func*(1d0/gamma**2+d2phi(3)))*wght
               Cmat2(n,D2_GAA_GAA) = Cmat2(n,D2_GAA_GAA) + t
               Cmat2(n,D2_GAA_GBB) = Cmat2(n,D2_GAA_GBB) + t
               Cmat2(n,D2_GAA_GAB) = Cmat2(n,D2_GAA_GAB) + t*2d0
               Cmat2(n,D2_GAB_GAB) = Cmat2(n,D2_GAB_GAB) + t*4d0
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
#ifdef SECOND_DERIV
            d2rs = -FOUR3*d1rs*arho
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
#ifdef SECOND_DERIV
            rrho2 = 2.d0*arho*arho
c           1 = aa, 2 = ab, 3 = bb
            d2z(1) =-rrho2*(1.d0-zeta)
            d2z(2) = rrho2*zeta
            d2z(3) = rrho2*(1.d0+zeta)
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
#ifdef SECOND_DERIV
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
c         
c           C(n)
c         
            anum = fff+alpha*rs+beta*rs2
            aden = 1.d0+pgamma*rs+delta*rs2+beta10*rs3
            Cn = zzz + anum/aden
            d1anum = alpha + 2d0*beta*rs
            d1aden = pgamma + 2d0*delta*rs + 3d0*beta10*rs2
#ifdef SECOND_DERIV
            d2anum = 2d0*beta
            d2aden = 2d0*delta + 6d0*beta10*rs
#endif
c     First compute rs derivative
            d1Cn = d1anum/aden - anum*d1aden/aden**2
#ifdef SECOND_DERIV
            d2Cn = d2anum/aden - (2d0*d1anum*d1aden+anum*d2aden)/aden**2
     &           + 2d0*anum*d1aden**2/aden**3
#endif
c     Convert to rho derivative
#ifdef SECOND_DERIV
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
#ifdef SECOND_DERIV
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
#ifdef SECOND_DERIV
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
#include "nwpwxc_c_perdew86.F"
#endif
C>
C> @}

